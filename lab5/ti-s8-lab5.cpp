#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>

using namespace std;

void encode_Hamming(string filename)
{
    ifstream input(filename);
    ofstream output("encoded_Hamming_" + filename);

    char c;
    vector <int> buffer;
    while (input.get(c)) {
        buffer.push_back(c - '0');

        if (buffer.size() == 4) {
            int encoded_buffer[7]{};
            encoded_buffer[2] = buffer[0];
            encoded_buffer[4] = buffer[1];
            encoded_buffer[5] = buffer[2];
            encoded_buffer[6] = buffer[3];

            encoded_buffer[0] = encoded_buffer[2] ^ encoded_buffer[4] ^ encoded_buffer[6];
            encoded_buffer[1] = encoded_buffer[2] ^ encoded_buffer[5] ^ encoded_buffer[6];
            encoded_buffer[3] = encoded_buffer[4] ^ encoded_buffer[5] ^ encoded_buffer[6];

            for (int i = 0; i < 7; i++)
                output.put(encoded_buffer[i] + '0');

            buffer.clear();
        }
    }

    input.close();
    output.close();
}

void change_bits_with_probability(string filename, double probability)
{
    ifstream input(filename);
    ofstream output("probability_" + to_string(probability) + "_" + filename);
    
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> distribution(0.0, 1.0);

    char c;
    while (input.get(c)) {
        if (distribution(gen) <= probability)
            c = c == '0' ? '1' : '0';
        output.put(c);
    }

    input.close();
    output.close();
}

void decode_Hamming(string filename)
{
    ifstream input(filename);
    ofstream output("decoded_Hamming_" + filename);

    long error_count = 0;

    char c;
    vector <int> buffer;
    while (input.get(c)) {
        buffer.push_back(c - '0');

        if (buffer.size() == 7) {
            int data[4]{};

            int p1 = buffer[0] ^ buffer[2] ^ buffer[4] ^ buffer[6];
            int p2 = buffer[1] ^ buffer[2] ^ buffer[5] ^ buffer[6];
            int p3 = buffer[3] ^ buffer[4] ^ buffer[5] ^ buffer[6];

            int error_bit = p1 + 2 * p2 + 4 * p3;
            if (error_bit != 0) {
                buffer[error_bit - 1] = buffer[error_bit - 1] == 0 ? 1 : 0;
                error_count++;
            }

            data[0] = buffer[2];
            data[1] = buffer[4];
            data[2] = buffer[5];
            data[3] = buffer[6];

            for (int i = 0; i < 4; i++)
                output.put(data[i] + '0');

            buffer.clear();
        }
    }

    cout << "Error bits found during decoding = " << error_count << "\n";

    input.close();
    output.close();
}

void compare_original_and_decoded_files(string original_filename, string decoded_filename)
{
    ifstream original(original_filename);
    ifstream decoded(decoded_filename);

    long error_count = 0;

    char original_char, decoded_char;
    while (original.get(original_char) && decoded.get(decoded_char))
        if (original_char != decoded_char)
            error_count++;

    cout << "Errors found when comparing files = " << error_count << "\n";

    original.close();
    decoded.close();
}

int main()
{
    string filename = "encoded_Huffman_Beauty_and_the_Beast.txt";
    // orignal file = 149248 bits
    string encoded_filename = "encoded_Hamming_" + filename;

    encode_Hamming(filename);

    vector<double> probabilities = { 0.0001, 0.001, 0.01, 0.1 };
    
    for (double probability : probabilities) {
        cout << "Probability = " << probability << "\n";
        change_bits_with_probability(encoded_filename, probability);

        string probability_encoded_filename = "probability_" + to_string(probability) + "_" + encoded_filename;
        decode_Hamming(probability_encoded_filename);

        string decoded_probability_encoded_filename = "decoded_Hamming_" + probability_encoded_filename;
        compare_original_and_decoded_files(filename, decoded_probability_encoded_filename);
        cout << "\n";
    }
    
    return 0;
}
