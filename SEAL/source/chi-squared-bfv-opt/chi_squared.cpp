#include "chi_squared.h"

#include "../common.h"

void ChiSquared::setup_context_bfv(std::size_t poly_modulus_degree,
                                   std::uint64_t plain_modulus) {
    /// Wrapper for parameters
    seal::EncryptionParameters params(seal::scheme_type::bfv);
    params.set_poly_modulus_degree(poly_modulus_degree);
    // Cingulata params + an additional moduli (44) as computation otherwise
    // cannot be performed

#ifdef MANUALPARAMS
    params.set_coeff_modulus(seal::CoeffModulus::Create(
            poly_modulus_degree, {60, 60, 30}));
#endif

#ifdef CINGUPARAM
    params.set_coeff_modulus(seal::CoeffModulus::Create(
        poly_modulus_degree, {30, 40, 44, 50, 54, 60, 60}));
#endif

#ifdef SEALPARAMS
    params.set_coeff_modulus(seal::CoeffModulus::BFVDefault(
        poly_modulus_degree, seal::sec_level_type::tc128));
#endif
    // set plaintext modulus suitable for batching
    params.set_plain_modulus(seal::PlainModulus::Batching(poly_modulus_degree, 20));
    // Instantiate context
    context = std::make_shared<seal::SEALContext>(params);

    /// Create keys
    seal::KeyGenerator keyGenerator(*context);
    keyGenerator.create_public_key(publicKey);
    secretKey = keyGenerator.secret_key();
    keyGenerator.create_relin_keys(relinKeys);

    // Provide both public and secret key, however, we will use public-key
    // encryption as this is the one used in a typical client-server scenario.
    encryptor =
            std::make_unique<seal::Encryptor>(*context, publicKey, secretKey);
    evaluator = std::make_unique<seal::Evaluator>(*context);
    decryptor = std::make_unique<seal::Decryptor>(*context, secretKey);
    encoder = std::make_unique<seal::BatchEncoder>(*context);

    auto qualifiers = context->first_context_data()->qualifiers();
    std::cout << "Batching enabled: " << std::boolalpha << qualifiers.using_batching << std::endl;
}

namespace {
    void log_time(std::stringstream &ss,
                  std::chrono::time_point<std::chrono::high_resolution_clock> start,
                  std::chrono::time_point<std::chrono::high_resolution_clock> end,
                  bool last = false) {
        ss << std::chrono::duration_cast<ms>(end - start).count();
        if (!last) ss << ",";
    }
}  // namespace

uint64_t ChiSquared::get_decrypted_value(seal::Ciphertext value) {
    seal::Plaintext tmp;
    std::vector<uint64_t> resultvec(encoder->slot_count(), 0ULL);;
    decryptor->decrypt(value, tmp);
    encoder->decode(tmp, resultvec);
    return resultvec[0];
}

ResultCiphertexts ChiSquared::compute_alpha_betas(const seal::Ciphertext &N_0,
                                                  const seal::Ciphertext &N_1,
                                                  const seal::Ciphertext &N_2) {
    std::size_t slot_count = encoder->slot_count();
    std::vector<uint64_t> four_vec(slot_count, 4ULL);
    seal::Plaintext four;
    encoder->encode(four_vec, four);

    std::vector<uint64_t> two_vec(slot_count, 4ULL);
    seal::Plaintext two;
    encoder->encode(two_vec, two);

    // compute alpha
    std::cout << "Computing alpha" << std::endl;
    seal::Ciphertext alpha;
    evaluator->multiply_plain(N_0, four, alpha);
    evaluator->relinearize_inplace(alpha, relinKeys);
    evaluator->multiply_inplace(alpha, N_2);
    evaluator->relinearize_inplace(alpha, relinKeys);
    seal::Ciphertext N_1_pow2;
    evaluator->exponentiate(N_1, 2, relinKeys, N_1_pow2);
    evaluator->sub_inplace(alpha, N_1_pow2);
    evaluator->exponentiate_inplace(alpha, 2, relinKeys);

    // compute beta_1
    seal::Ciphertext beta_1;
    seal::Ciphertext N_0_t2;
    evaluator->multiply_plain(N_0, two, N_0_t2);
    evaluator->relinearize_inplace(N_0_t2, relinKeys);
    seal::Ciphertext twot_N_0__plus__N_1;
    evaluator->add(N_0_t2, N_1, twot_N_0__plus__N_1);
    evaluator->exponentiate(twot_N_0__plus__N_1, 2, relinKeys, beta_1);
    evaluator->multiply_plain_inplace(beta_1, two);
    evaluator->relinearize_inplace(beta_1, relinKeys);

    // compute beta_2
    seal::Ciphertext beta_2;
    seal::Ciphertext t2_N_2;
    evaluator->multiply_plain(N_2, two, t2_N_2);
    evaluator->relinearize_inplace(t2_N_2, relinKeys);
    seal::Ciphertext twot_N_2__plus__N_1;
    evaluator->add(t2_N_2, N_1, twot_N_2__plus__N_1);
    evaluator->multiply(twot_N_0__plus__N_1, twot_N_2__plus__N_1, beta_2);
    evaluator->relinearize_inplace(beta_2, relinKeys);

    // compute beta_3
    seal::Ciphertext beta_3;
    evaluator->exponentiate(twot_N_2__plus__N_1, 2, relinKeys, beta_3);
    evaluator->multiply_plain_inplace(beta_3, two);
    evaluator->relinearize_inplace(beta_3, relinKeys);

    return ResultCiphertexts(alpha, beta_1, beta_2, beta_3);
}

void ChiSquared::run_chi_squared() {
    std::stringstream ss_time;

    // set up the BFV scheme
    auto t0 = Time::now();
    setup_context_bfv(32768, 4096);
    auto t1 = Time::now();
    log_time(ss_time, t0, t1, false);

    auto t2 = Time::now();
    int32_t n0_val = 2, n1_val = 7, n2_val = 9;
    seal::Ciphertext n0, n1, n2;
    seal::Plaintext n0_plain;
    seal::Plaintext n1_plain;
    seal::Plaintext n2_plain;

    //initialise vectors for each integer (we use batch encoding)
    std::size_t slot_count = encoder->slot_count();
    std::vector<uint64_t> n0_vec(slot_count, n0_val);
    std::vector<uint64_t> n1_vec(slot_count, n1_val);
    std::vector<uint64_t> n2_vec(slot_count, n2_val);
    encoder->encode(n0_vec, n0_plain);
    encoder->encode(n1_vec, n1_plain);
    encoder->encode(n2_vec, n2_plain);
    encryptor->encrypt(n0_plain, n0);
    encryptor->encrypt(n1_plain, n1);
    encryptor->encrypt(n2_plain, n2);
    auto t3 = Time::now();
    log_time(ss_time, t2, t3, false);

    // perform FHE computation
    auto t4 = Time::now();
    auto result = compute_alpha_betas(n0, n1, n2);
    auto t5 = Time::now();
    log_time(ss_time, t4, t5, false);

    // decrypt results
    auto t6 = Time::now();
    int32_t result_alpha = get_decrypted_value(result.alpha);
    int32_t result_beta1 = get_decrypted_value(result.beta_1);
    int32_t result_beta2 = get_decrypted_value(result.beta_2);
    int32_t result_beta3 = get_decrypted_value(result.beta_3);
    auto t7 = Time::now();
    log_time(ss_time, t6, t7, true);

    // check results
    auto exp_alpha = std::pow((4 * n0_val * n2_val) - std::pow(n1_val, 2), 2);
    assert(("Unexpected result for 'alpha' encountered!",
            result_alpha == exp_alpha));
    std::cout << "Expected alpha: " << exp_alpha << ", calculated alpha: " << result_alpha << std::endl;
    auto exp_beta_1 = 2 * std::pow(2 * n0_val + n1_val, 2);
    assert(("Unexpected result for 'beta_1' encountered!",
            result_beta1 == exp_beta_1));
    std::cout << "Expected beta_1: " << exp_beta_1 << ", calculated beta_1: " << result_beta1 << std::endl;
    auto exp_beta_2 = ((2 * n0_val) + n1_val) * ((2 * n2_val) + n1_val);
    assert(("Unexpected result for 'beta_2' encountered!",
            result_beta2 == exp_beta_2));
    std::cout << "Expected beta_2: " << exp_beta_2 << ", calculated beta_2: " << result_beta2 << std::endl;
    auto exp_beta_3 = 2 * std::pow(2 * n2_val + n1_val, 2);
    assert(("Unexpected result for 'beta_3' encountered!",
            result_beta3 == exp_beta_3));
    std::cout << "Expected beta_3: " << exp_beta_3 << ", calculated beta_3: " << result_beta3 << std::endl;
    // write ss_time into file
    std::ofstream myfile;
    auto out_filename = std::getenv("OUTPUT_FILENAME");
    myfile.open(out_filename, std::ios_base::app);
    myfile << ss_time.str() << std::endl;
    myfile.close();

    // write FHE parameters into file
    write_parameters_to_file(context, "fhe_parameters_chi_squared.txt");
}

int main(int argc, char *argv[]) {
    std::cout << "Starting benchmark 'chi-squared-bfv-opt'..." << std::endl;
    ChiSquared().run_chi_squared();
    return 0;
}
