#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <vector>
#include <utility>
#include <queue>

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

void check_probabilities(vector<pair<char, double>> probabilities) {
    double remains = 1;
    int index_of_max = 0;

    for (int i = 0; i < probabilities.size(); i++) {
        remains -= probabilities[i].second;
        if (probabilities[i].second > probabilities[index_of_max].second)
            index_of_max = i;
    }

    probabilities[index_of_max].second += remains;
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

vector<pair<char, double>> calculate_probabilities(string file_name) {
    map<char, long> alphabet;
    long counter = 0;

    ifstream file(file_name);

    char c;
    while (file.get(c)) {
        if (check_char(c) != -1) {
            alphabet[check_char(c)]++;
            counter++;
        }
    }

    file.close();

    vector<pair<char, double>> probabilities;

    for (const auto& symbol : alphabet) {
        probabilities.push_back(make_pair(symbol.first, (double)symbol.second / counter));
        //cout << symbol.first << " " << symbol.second << "\n";
    }

    check_probabilities(probabilities);

    return probabilities;
}

double calculate_average_code_length(map<char, string> codes, map<char, double> probabilities) {
    double length = 0;
    for (const auto& code : codes)
        length += code.second.size() * probabilities.at(code.first);
    return length;
}

void encode_file(string file_name, string encoding_name, map<char, string> codes, map<char, double> probabilities) {
    ifstream input(file_name);
    ofstream output("encoded_" + encoding_name + "_" + file_name);

    char c;
    while (input.get(c)) {
        if (check_char(c) != -1) {
            output << codes.at(check_char(c));
        }
    }

    input.close();
    output.close();
}

void evaluate_file_encoding(string file_name, string encoding_name, map<char, string> codes, map<char, double> probabilities) {
    double avg_length = calculate_average_code_length(codes, probabilities);

    double max_entropy = -1;
    double entropy[3];
    for (int i = 0; i < 3; i++) {
        entropy[i] = get_evaluation("encoded_" + encoding_name + "_" + file_name, i + 1);
        max_entropy = max(max_entropy, entropy[i]);
    }

    cout << "FILE : " << file_name << "\n";
    cout << "ENCODING NAME : " << encoding_name << "\n";
    cout << "\tAVERAGE CODE LENGTH : " << avg_length << "\n";
    
    for (int i = 0; i < 3; i++) {
        cout << "\tENTROPY OF " << i + 1 << " : " << entropy[i] << "\n";
        //cout << "\t\tCODING REDUNDANCY : " << avg_length - entropy[i] << "\n";
    }

    cout << "\tCODING REDUNDANCY : " << avg_length - max_entropy << "\n\n";
    
    /*
    for (const auto& code : codes)
        cout << code.first << " " << code.second << " | " << probabilities[code.first] << "\n";   
    cout << "\n";
    */
}

struct Huffman_node {
    char symbol;
    double probability;
    Huffman_node* left;
    Huffman_node* right;

    Huffman_node(char symbol, double probability) : symbol(symbol), probability(probability), left(nullptr), right(nullptr) {}
};

struct Compare_nodes {
    bool operator()(Huffman_node* lhs, Huffman_node* rhs) const {
        return lhs->probability > rhs->probability;
    }
};

Huffman_node* build_Huffman_tree(vector<pair<char, double>>& probabilities) {
    priority_queue<Huffman_node*, vector<Huffman_node*>, Compare_nodes> node_queue;

    for (const auto& probability : probabilities) {
        node_queue.push(new Huffman_node(probability.first, probability.second));
    }

    while (node_queue.size() > 1) {
        Huffman_node* left_node = node_queue.top();
        node_queue.pop();

        Huffman_node* right_node = node_queue.top();
        node_queue.pop();

        Huffman_node* merged_node = new Huffman_node('\0', left_node->probability + right_node->probability);
        merged_node->left = left_node;
        merged_node->right = right_node;

        node_queue.push(merged_node);
    }

    return node_queue.top();
}

void build_Huffman_codes(Huffman_node* root, const string& code, map<char, string>& huffman_codes) {
    if (root) {
        if (!root->left && !root->right) {
            huffman_codes[root->symbol] = code;
            return;
        }

        build_Huffman_codes(root->left, code + "0", huffman_codes);
        build_Huffman_codes(root->right, code + "1", huffman_codes);
    }
}

void encode_Huffman_file(string file_name) {
    vector<pair<char, double>> probabilities = calculate_probabilities(file_name);
    map<char, string> huffman_codes;

    Huffman_node* root = build_Huffman_tree(probabilities);
    build_Huffman_codes(root, "", huffman_codes);

    map<char, double> probabilities_map;
    for (const auto& probability : probabilities)
        probabilities_map[probability.first] = probability.second;

    encode_file(file_name, "Huffman", huffman_codes, probabilities_map);
    evaluate_file_encoding(file_name, "Huffman", huffman_codes, probabilities_map);
}

void build_Fano_codes(vector<pair<char, double>>& probabilities, int start, int end, string code, map<char, string>& fano_codes) {
    if (start == end) {
        fano_codes[probabilities[start].first] = code;
        return;
    }

    int mid = start;
    double min_delta = 10;
    double left_sum = 0;
    double right_sum = 0;

    for (int i = start; i <= end; i++)
        right_sum += probabilities[i].second;

    for (int i = start; i <= end; i++) {
        left_sum += probabilities[i].second;
        right_sum -= probabilities[i].second;

        if (abs(left_sum - right_sum) < min_delta) {
            mid = i;
            min_delta = abs(left_sum - right_sum);
        }
    }

    build_Fano_codes(probabilities, start, mid, code + '0', fano_codes);
    build_Fano_codes(probabilities, mid + 1, end, code + '1', fano_codes);
}

void encode_Fano_file(string file_name) {
    vector<pair<char, double>> probabilities = calculate_probabilities(file_name);
    sort(probabilities.begin(), probabilities.end(),
        [](const auto& left, const auto& right) { return left.second > right.second; });
    
    map<char, string> fano_codes;
    build_Fano_codes(probabilities, 0, probabilities.size() - 1, "", fano_codes);
    
    map<char, double> probabilities_map;
    for (const auto& probability : probabilities)
        probabilities_map[probability.first] = probability.second;
    
    encode_file(file_name, "Fano", fano_codes, probabilities_map);
    evaluate_file_encoding(file_name, "Fano", fano_codes, probabilities_map);
}

int main() {
    vector<char> alphabet = { 'a', 'b', 'c', 'd', 'e', 'f'};
    vector<double> alphabet_probabilities = { 0.36, 0.181, 0.179, 0.12, 0.09, 0.07 };

    if (alphabet.size() != alphabet_probabilities.size()) {
        cout << "ERROR";
        return -1;
    }
    check_probabilities(alphabet_probabilities);

    generate_file_equal_prob(alphabet);
    generate_file_diff_prob(alphabet, alphabet_probabilities);
    
    vector<string> files = { "equal_prob.txt" , "diff_prob.txt" , "Beauty_and_the_Beast.txt" };
    for (string file : files) {
        encode_Fano_file(file);
        encode_Huffman_file(file);
    }

    return 0;
}