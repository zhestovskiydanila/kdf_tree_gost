#include "kdf_gost.hpp"

HMAC256::HMAC256() : key{} {};

HMAC256::HMAC256(std::vector<uint8_t> key_v) : key{} {
  if (key_v.size() < 32 || key_v.size() > 64) {
    perror("Неверная длина ключа для HMAC\n");
    exit(EXIT_FAILURE);
  }
  std::copy(key_v.begin(), key_v.end(), key.begin());
}

HMAC256::~HMAC256() { flush_key(key); }

void HMAC256::flush_key(hashb_t &key) {
  for (int i = 0; i < BLOCK_SIZE; ++i) {
    key[i] = static_cast<uint8_t>(dist(rd));
  }
}

void HMAC256::XOR512(hashb_t &lhs, const hashb_t &rhs) {
  for (int i = 0; i < BLOCK_SIZE; ++i) {
    lhs[i] ^= rhs[i];
  }
}

std::random_device HMAC256::rd;
std::uniform_int_distribution<uint8_t> HMAC256::dist(0, 255);

hash256_t HMAC256::generate_HMAC(std::vector<uint8_t> data) {
  Streebog256 hash{};
  hashb_t key_ipad{this->key}, key_opad{this->key};
  std::cout << "Key\n" << this->key << "\n";
  XOR512(key_ipad, ipad);
  XOR512(key_opad, opad);
  std::cout << "Key_ipad\n" << key_ipad << "\nData\n";
  for (int i = 0; i < data.size(); ++i) {
    std::cout << std::setfill('0') << std::setw(2) << std::hex << static_cast<unsigned int>(data[i]);
  }
  std::cout << "\n";
  std::vector<uint8_t> first_hash_data{};
  first_hash_data.reserve(BLOCK_SIZE + data.size());
  first_hash_data.insert(first_hash_data.end(), key_ipad.begin(),
                         key_ipad.end());
  first_hash_data.insert(first_hash_data.end(), data.begin(), data.end());
  std::cout << "Concat key and data\n";
  for (int i = 0; i < first_hash_data.size(); ++i) {
    std::cout << std::setfill('0') << std::setw(2) << std::hex << static_cast<unsigned int>(first_hash_data[i]);
  }
  std::cout << "\n";
  hash256_t first_hash =
      hash.digest_message(first_hash_data.data(), first_hash_data.size());
  std::cout << "First hash\n" << first_hash << "\n";
  std::vector<uint8_t> second_hash_data{};
  second_hash_data.reserve(BLOCK_SIZE + HASH_SIZE);
  second_hash_data.insert(second_hash_data.end(), key_opad.begin(),
                          key_opad.end());
  second_hash_data.insert(second_hash_data.end(), first_hash.begin(),
                          first_hash.end());
  return hash.digest_message(second_hash_data.data(), second_hash_data.size());
}

int main() {
  std::vector<uint8_t> key = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
                              0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
                              0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
                              0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};
  std::vector<uint8_t> data = {0x01, 0x26, 0xbd, 0xb8, 0x78, 0x00, 0xaf, 0x21,
                               0x43, 0x41, 0x45, 0x65, 0x63, 0x78, 0x01, 0x00};

  hash256_t result = {0xa1, 0xaa, 0x5f, 0x7d, 0xe4, 0x02, 0xd7, 0xb3,
                      0xd3, 0x23, 0xf2, 0x99, 0x1c, 0x8d, 0x45, 0x34,
                      0x01, 0x31, 0x37, 0x01, 0x0a, 0x83, 0x75, 0x4f,
                      0xd0, 0xaf, 0x6d, 0x7c, 0xd4, 0x92, 0x2e, 0xd9};

  HMAC256 hmac(key);
  hash256_t hmac_result = hmac.generate_HMAC(data);
  if (result != hmac_result) {
    std::cout << "True HMAC\n" << result << "\n" << "Test HMAC\n" << hmac_result << "\n";
  }
  return 0;
}