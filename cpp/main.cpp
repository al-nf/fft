#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
using namespace std;

int log2(int N) {
    int i = 0;
    while (N > 1) {
        N >>= 1;
        i++;
    }
    return i;
}

int reverse(int N, int n) {
    int result = 0;
    int bits = log2(N);
    for (int i = 0; i < bits; i++) {
        if (n & (1 << i)) {
            result |= (1 << (bits - 1 - i));
        }
    }
    return result;
}

void transform(vector<complex<double> >& f) {
    int N = f.size();

    for (int i = 0; i < N; i++) {
        int j = reverse(N, i);
        if (i < j) {
            swap(f[i], f[j]);
        }
    }

    vector<complex<double> > W(N / 2);
    for (int i = 0; i < N / 2; i++) {
        W[i] = polar(1.0, -2.0 * M_PI * i / N);
    }

    // FFT
    for (int len = 2; len <= N; len *= 2) {
        int halfLen = len / 2;
        for (int i = 0; i < N; i += len) {
            for (int j = 0; j < halfLen; j++) {
                complex<double> u = f[i + j];
                complex<double> v = W[N / len * j] * f[i + j + halfLen];
                f[i + j] = u + v;
                f[i + j + halfLen] = u - v;
            }
        }
    }
}

void FFT(vector<complex<double> >& f, double d) {
    transform(f);
    for (auto& val : f) {
        val *= d;
    }
}

int main() {
    vector<complex<double> > f;
    f.push_back(complex<double>(1, 0));
    f.push_back(complex<double>(2, 0));
    f.push_back(complex<double>(3, 0));
    f.push_back(complex<double>(4, 0));
    f.push_back(complex<double>(0, 0));
    f.push_back(complex<double>(0, 0));
    f.push_back(complex<double>(0, 0));
    f.push_back(complex<double>(0, 0));

    FFT(f, 1.0);

    for (const auto& val : f) {
        cout << val << '\n';
    }
}
