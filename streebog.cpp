#include "kdf_gost.hpp"
#include <vector>

hashb_t Streebog256::func_E(hashb_t hash, hashb_t msg, hashb_t N) {
  std::array<hashb_t, 13> Keys{};
  hashb_t E_hash = msg;
  XOR512(hash, N);
  LPS(hash);
  Keys[0] = hash;
  for (int i = 0; i < 12; ++i) {
    hashb_t tmp = Keys[i];
    XOR512(tmp, C[i]);
    LPS(tmp);
    Keys[i + 1] = tmp;
  }
  for (int i = 0; i < 12; ++i) {
    XOR512(E_hash, Keys[i]);
    LPS(E_hash);
  }

  XOR512(E_hash, Keys[12]);

  return E_hash;
}

void Streebog256::XOR512(hashb_t &lhs, const hashb_t &rhs) {
#ifdef __AVX512F__
  __m512i acc = _mm512_setzero_si512();
  acc = _mm512_xor_si512(*reinterpret_cast<__m512i *>(lhs.data()),
                         *reinterpret_cast<const __m512i *>(rhs.data()));
  _mm512_store_si512(reinterpret_cast<__m512 *>(lhs.data()), acc);
#else
  uint64_t *vec_lhs = reinterpret_cast<uint64_t *>(lhs.data());
  const uint64_t *vec_rhs = reinterpret_cast<const uint64_t *>(rhs.data());
  for (int i = 0; i < 8; ++i) {
    vec_lhs[i] ^= vec_rhs[i];
  }
#endif
}

void Streebog256::ADD512(hashb_t &lhs, const hashb_t &rhs) {
  int carry = 0;
  for (int i = BLOCK_SIZE - 1; i >= 0; --i) {
    int sum = lhs[i] + rhs[i] + carry;
    lhs[i] = sum & 0xFF;
    carry = sum >> 8;
  }
}

void Streebog256::ADD512_uint(hashb_t &lhs, uint64_t rhs) {
  hashb_t hash_rhs{};
  uint64_t rhs_le = __builtin_bswap64(rhs);
  uint64_t *vec_rhs = reinterpret_cast<uint64_t *>(hash_rhs.data());
  vec_rhs[7] = rhs_le;

  ADD512(lhs, hash_rhs);
}

void Streebog256::func_S(hashb_t &msg) {
#ifdef __AVX512BW__
  __m512i acc = _mm512_load_epi64(reinterpret_cast<__m512i *>(msg.data()));
  __m512i res = _mm512_shuffle_epi8(Pi_acc, acc);
  _mm512_storeu_epi64(reinterpret_cast<__m512i *>(msg.data()), res);
#else
  for (int i = 0; i < BLOCK_SIZE; ++i) {
    msg[i] = Pi[msg[i]];
  }
#endif
}

void Streebog256::func_P(hashb_t &msg) {
#ifdef __AVX512F__
  __m512i p_reg =
      _mm512_permutexvar_epi8(revTau, *reinterpret_cast<__m512i *>(msg.data()));
  _mm512_storeu_epi8(reinterpret_cast<__m512i *>(msg.data()), p_reg);
#else
  hashb_t tmp{};
  for (int i = 0; i < BLOCK_SIZE; ++i) {
    tmp[i] = msg[Tau[i]];
  }
  msg = tmp;
#endif
}

void Streebog256::func_L(hashb_t &msg) {
#ifdef __AVX512F__
  __m512i acc = _mm512_setzero_si512();
  __m512i shift_acc =
      _mm512_loadu_si512(reinterpret_cast<__m512 *>(msg.data()));
  for (int i = 0; i < BLOCK_SIZE; ++i) {
    shift_acc = _mm512_srli_epi64(shift_acc, 1);
    __mmask8 mask = _mm512_movepi64_mask(shift_acc);
    __m512i xor_reg = _mm512_maskz_set1_epi64(mask, A_le[63 - i]);
    acc = _mm512_xor_epi64(acc, xor_reg);
  }
  _mm512_storeu_si512(reinterpret_cast<__m512i *>(msg.data()), acc);
#else
  hashb_t tmp{};
  uint64_t *tmp_vec = reinterpret_cast<uint64_t *>(tmp.data());
  const uint64_t *vec = reinterpret_cast<uint64_t *>(msg.data());
  for (int i = 0; i < 8; ++i) {
    uint64_t num = __builtin_bswap64(vec[i]);
    for (int j = 0; j < 64; ++j) {
      if ((num >> j) & 1) {
        tmp_vec[i] ^= A_le[63 - j];
      }
    }
  }
  msg = tmp;
#endif
}

void Streebog256::LPS(hashb_t &msg) {
  func_S(msg);

  func_P(msg);

  func_L(msg);
}

hashb_t Streebog256::pad_block(const hashb_t &msg, size_t msg_len) {
  hashb_t padded_block{};
  int offset = BLOCK_SIZE - msg_len;
  if (offset > 0) {
    memcpy(padded_block.data(), msg.data(), BLOCK_SIZE);
    padded_block[offset - 1] = 0x1;
  } else {
    padded_block = msg;
  }

  return padded_block;
}

void Streebog256::func_gN(hashb_t &hash, hashb_t msg, hashb_t N) {
  hashb_t E_hash = hash;

  E_hash = func_E(E_hash, msg, N);

  XOR512(hash, msg);
  XOR512(hash, E_hash);
}

void Streebog256::func_g0(hashb_t &hash, hashb_t msg) {
  hashb_t E_hash = hash;
  hashb_t N{};
  E_hash = func_E(E_hash, msg, N);
  XOR512(hash, msg);
  XOR512(hash, E_hash);
}

hash256_t Streebog256::digest_message(uint8_t *msg, int64_t len) {
  hashb_t hash{init_vector};
  hashb_t N{};
  hashb_t Sigma{};
  int64_t pos = len;
  while (len >= BLOCK_SIZE) {
    pos -= BLOCK_SIZE;
    hashb_t hash_msg{};
    memcpy(hash_msg.data(), msg + pos, BLOCK_SIZE);
    len -= BLOCK_SIZE;
    func_gN(hash, hash_msg, N);
    ADD512_uint(N, 512);

    ADD512(Sigma, hash_msg);
  }
  hashb_t hash_msg{};
  if (len > 0) {
    memcpy(hash_msg.data() + BLOCK_SIZE - len, msg, len);
    hash_msg = pad_block(hash_msg, len);
  } else {
    hash_msg = pad_block(hash_msg, 0);
  }
  func_gN(hash, hash_msg, N);
  if (len > 0) {
    ADD512_uint(N, static_cast<uint64_t>(len * 8));
  } else {
    ADD512_uint(N, static_cast<uint64_t>(BLOCK_SIZE * 8));
  }

  ADD512(Sigma, hash_msg);

  func_g0(hash, N);
  func_g0(hash, Sigma);
  hash256_t final_hash{};
  memcpy(final_hash.data(), hash.data(), HASH_SIZE);
  return final_hash;
}

std::vector<uint8_t> hexToBytes(std::string &hex) {

  std::vector<uint8_t> bytes;
  bytes.reserve(hex.size() / 2);

  for (size_t i = 0; i < hex.size(); i += 2) {
    std::string byteStr = hex.substr(i, 2);
    uint8_t byte = static_cast<uint8_t>(std::stoul(byteStr, nullptr, 16));
    bytes.push_back(byte);
  }

  return bytes;
}

std::array<uint8_t, 64> blockFromHex(std::string &hex) {
  auto vec = hexToBytes(hex);
  std::array<uint8_t, 64> arr{};
  std::copy(vec.begin(), vec.end(), arr.begin());
  return arr;
}

/* int main() {
  Streebog256 hash;

  std::string str2 =
      "fbe2e5f0eee3c820fbeafaebef20fffbf0e1e0f0f520e0ed20e8ece0eb"
      "e5f0f2f120fff0eeec20f120faf2fee5e2202ce8f6f3ede220e8e6eee1"
      "e8f0f2d1202ce8f0f2e5e220e5d1";

  std::string str1 =
      "32313039383736353433323130393837363534333231303938373635343332"
      "3130393837363534333231303938373635343332313039383736353433323130";

  std::vector<uint8_t> bytes = hexToBytes(str1);

  std::string true_hash1 =
      "00557be5e584fd52a449b16b0251d05d27f94ab76cbaa6da890b59d8ef1e159d";

  std::vector<uint8_t> true_hashv1 = hexToBytes(true_hash1);

  std::string true_hash2 =
      "508f7e553c06501d749a66fc28c6cac0b005746d97537fa85d9e40904efed29d";

  std::vector<uint8_t> true_hashv2 = hexToBytes(true_hash2);
  std::array<uint8_t, 32> true_hasharr1{};
  std::array<uint8_t, 32> true_hasharr2{};
  std::copy(true_hashv1.begin(), true_hashv1.end(), true_hasharr1.begin());
  std::copy(true_hashv2.begin(), true_hashv2.end(), true_hasharr2.begin());

  std::vector<uint8_t> bytes2 = hexToBytes(str2);
  hash256_t test_hash1 = hash.digest_message(bytes.data(), bytes.size());

  std::cout << "Test hash\n"
            << test_hash1 << "\n"
            << "True hash\n"
            << true_hash1 << "\n";

  hash256_t test_hash2 = hash.digest_message(bytes2.data(), bytes2.size());

  std::cout << "Test hash\n"
            << test_hash2 << "\n"
            << "True hash\n"
            << true_hash2 << "\n";
  std::string test_str =
      "fbeafaebef20fffbf0e1e0f0f520e0ed20e8ece0ebe5f0f2f120fff0eeec20f120faf2fe"
      "e5e2202ce8f6f3ede220e8e6eee1e8f0f2d1202ce8f0f2e5e220e5d1";
  std::vector<uint8_t> test_vec = hexToBytes(test_str);
  std::array<uint8_t, 64> test_array{};
  std::copy(test_vec.begin(), test_vec.end(), test_array.begin());

  std::string test_str2 =
      "8d49118311e4d9e44fe2012b1faee26a9304dd7714cd311482ada7ad959fad00"
      "87c8475d0c0e2c0e47470abce8473847a73b4157572f57a56cd15b2d0bd20b86";

  std::string str2_after_S =
      "a3a72a2e0fb5e6f812681222fec037b0db972086a395a387a6084508cae13093"
      "aa71d352dcbce288e9a39718a727f6fd4c5da5d0bc10fac3707ccd127fe45475";

  auto arr1 = blockFromHex(test_str2);
  auto arr2 = blockFromHex(str2_after_S);
  return 0;
} */