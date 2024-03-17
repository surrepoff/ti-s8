#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <vector>

using namespace std;

char check_char(char c) {
    if (!isalpha(c) && !isdigit(c) && c != ' ')
        return -1;

    if (isalpha(c))
        return tolower(c);

    return c;
}

void check_probabilities(vector<double> probabilities) {
    double remains = 1;
    int index_of_max = 0;

    for (int i = 0; i < probabilities.size(); i++) {
        remains -= probabilities[i];
        if (probabilities[i] > probabilities[index_of_max])
            index_of_max = i;
    }

    probabilities[index_of_max] += remains;
}

void generate_file_equal_prob(vector<char> alphabet) {
    ofstream file("equal_prob.txt");

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, alphabet.size() - 1);

    for (long i = 0; i < 1 << 15; i++) {
        file << alphabet[distribution(gen)];
    }

    file.close();
}

void generate_file_diff_prob(vector<char> alphabet, vector<double> alphabet_probabilities) {
    ofstream file("diff_prob.txt");

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> distribution(alphabet_probabilities.begin(), alphabet_probabilities.end());

    for (long i = 0; i < 1 << 15; i++) {
        file << alphabet[distribution(gen)];
    }

    file.close();
}

double calculate_entropy(vector<double> probabilities) {
    double entropy = 0;
    for (double probability : probabilities)
        entropy += probability * log2(probability);

    return -entropy;
}

double get_evaluation(string file_name, int limit) {
    if (limit < 1 || limit > 100) return -1;

    map<string, long> alphabet;
    list<char> buffer;
    long counter = 0;

    ifstream file(file_name);

    char c;
    while (file.get(c)) {
        if (check_char(c) != -1) {
            if (buffer.size() == limit) {
                buffer.pop_front();
            }

            buffer.push_back(check_char(c));

            if (buffer.size() == limit) {
                string str = "";
                for (char symbol : buffer)
                    str += symbol;

                alphabet[str]++;
                counter++;
            }
        }
    }

    file.close();

    vector<double> probabilities;

    for (const auto& symbol : alphabet) {
        probabilities.push_back((double)symbol.second / counter);
        //cout << symbol.first << " " << symbol.second << "\n";
    }

    check_probabilities(probabilities);

    return calculate_entropy(probabilities) / limit;
}

int main() {
    vector<char> alphabet = { 'a', 'b', 'c', 'd', 'e' };
    vector<double> alphabet_probabilities = { 0.80, 0.05, 0.05, 0.05, 0.05 };

    if (alphabet.size() != alphabet_probabilities.size()) {
        cout << "ERROR";
        return -1;
    }
    check_probabilities(alphabet_probabilities);

    generate_file_equal_prob(alphabet);
    generate_file_diff_prob(alphabet, alphabet_probabilities);

    for (int i = 1; i <= 3; i++) {
        cout << "Entropy of " << i << " consecutive symbols:\n";
        cout << "equal prob: " << get_evaluation("equal_prob.txt", i) << "\n";
        cout << "diff prob: " << get_evaluation("diff_prob.txt", i) << "\n";
        cout << "text: " << get_evaluation("Beauty_and_the_Beast.txt", i) << "\n\n";
    }

    vector<double> alphabet_prob_equal;
    for (char c : alphabet)
        alphabet_prob_equal.push_back((double)1 / alphabet.size());
    check_probabilities(alphabet_prob_equal);

    cout << "Theoretical entropy:\n";
    cout << "equal prob: " << calculate_entropy(alphabet_prob_equal) << "\n";
    cout << "diff prob: " << calculate_entropy(alphabet_probabilities) << "\n";

    return 0;
}