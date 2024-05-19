#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "AudioFile.h"
#include <vector>
#include <cmath>
#include <random>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct wave {
    std::vector<double>x;
    std::vector<double>y;
};

wave sinus(double frequency) {
    wave sinus;
    for (int i = 0; i < 628; i++)
    {
        sinus.x.push_back(static_cast<double>(i) / 314);
        sinus.y.push_back(std::sin(3.14 * frequency * sinus.x[i]));
    }
    return sinus;
}
wave cosinus(double frequency) {
    wave cosinus;
    for (int i = 0; i < 628; i++)
    {
        cosinus.x.push_back(static_cast<double>(i) / 314);
        cosinus.y.push_back(std::cos(3.14 * frequency * cosinus.x[i]));
    }
    return cosinus;
}

wave sawtooth(int f, int fs, int num_periods) {
    wave sawtooth;
    double step = (2.0 * 32768) / (fs * f);
    double amplitude = 0;
    for (int p = 0; p < num_periods; ++p) {
        for (int i = 0; i < fs; ++i) {
            sawtooth.x.push_back(static_cast<double>(i + p * fs) / fs);
            sawtooth.y.push_back(amplitude);
            amplitude += step;
            if (amplitude >= 32768)
                amplitude -= 65536;
        }
    }
    return sawtooth;
}

wave square(int f, int fs, int num_periods) {
    wave square;
    for (int p = 0; p < num_periods; ++p) {
        for (int i = 0; i < fs; ++i) {
            square.x.push_back(static_cast<double>(i + p * fs) / fs);
            square.y.push_back(i < fs / 2 ? -32768 : 32768);
        }
    }
    return square;
}

std::vector<std::complex<double>> DFT(const wave& input_wave) {
    int N = input_wave.y.size();
    std::vector<std::complex<double>> X(N);
    for (int k = 0; k < N; ++k) {
        std::complex<double> sum = 0.0;
        for (int n = 0; n < N; ++n) {
            double angle = -2 * M_PI * k * n / N;
            sum += input_wave.y[n] * std::exp(std::complex<double>(0, angle));
        }
        X[k] = sum;
    }
    return X;
}

std::vector<double> IDFT(const std::vector<std::complex<double>>& X) {
    int N = X.size();
    std::vector<double> x(N);
    for (int n = 0; n < N; ++n) {
        std::complex<double> sum = 0.0;
        for (int k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            sum += X[k] * std::exp(std::complex<double>(0, angle));
        }
        x[n] = std::real(sum) / N;
    }

    return x;
}

void plot_idft(const std::vector<std::complex<double>>& X) {
    std::vector<double> time_domain_signal = IDFT(X);
    std::vector<double> x(time_domain_signal.size());
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = static_cast<double>(i);
    }
    matplot::plot(x, time_domain_signal);
    matplot::title("Time-Domain Signal (Inverse DFT)");
    matplot::xlabel("Sample Index");
    matplot::ylabel("Amplitude");
    matplot::show();
}

std::vector<double> magnitude(const std::vector<std::complex<double>>& vec) {
    std::vector<double> mag;
    for (const auto& val : vec) {
        mag.push_back(std::abs(val));
    }
    return mag;
}

void plot_dft(const wave& input_wave) {
    std::vector<std::complex<double>> dft_result = DFT(input_wave);
    std::vector<double> dft_magnitude = magnitude(dft_result);
    std::vector<double> dft_freq;
    int N = dft_result.size();
    for (int k = 0; k < N; ++k) {
        dft_freq.push_back(static_cast<double>(k) / N);
    }
    matplot::plot(dft_freq, dft_magnitude);
    matplot::title("Magnitude of DFT Coefficients");
    matplot::xlabel("Normalized Frequency");
    matplot::ylabel("Magnitude");
    matplot::show();
}

void plot_sinus( double frequency) {
    wave sine = sinus(frequency);
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ sinusoidalny");
    matplot::plot(sine.x, sine.y)->color({ 1.0f, 0.08f, 0.58f });
    matplot::show();
}

void plot_cosinus(double frequency) {
    wave cosine = cosinus(frequency);
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ cosinusoidalny");
    matplot::plot(cosine.x, cosine.y)->color({ 1.0f, 0.08f, 0.58f });;
    matplot::show();
}

void plot_sawtooth(int f, int fs, int num_periods) {
    wave saw = sawtooth(f, fs, num_periods);
    matplot::plot(saw.x, saw.y)->color({ 1.0f, 0.08f, 0.58f });
    matplot::xlabel("nr próbki");
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ pi³okszta³tny");
    matplot::show();
}

void plot_square(int f, int fs, int num_periods) {
    wave sq = square(f, fs, num_periods);
    matplot::plot(sq.x, sq.y)->color({ 1.0f, 0.08f, 0.58f });
    matplot::xlabel("nr próbki");
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ prostok¹tny");
    matplot::show();
}

void visualize_audio(const std::string& audioFilePath) {
    AudioFile<double> audioFile;
    bool load = audioFile.load(audioFilePath);
    if (!load) {
        std::cerr << "Audio is not loaded" << std::endl;
        return;
    }
    int samples = audioFile.getNumSamplesPerChannel();
    int channels = audioFile.getNumChannels();
    if (samples <= 0) {
        std::cerr << "audio has no samples" << std::endl;
        return;
    }
    std::vector<double> audioData(audioFile.samples[0].begin(), audioFile.samples[0].end());
    std::vector<double> x(samples);
    for (int i = 0; i < samples; ++i) {
        x[i] = static_cast<double>(i) / audioFile.getSampleRate();
    }
    matplot::plot(x, audioData)->color({ 1.0f, 0.08f, 0.58f });;
    matplot::title("Sygna³ Audio");
    matplot::show();
}

wave add_noise(const wave& input_wave, double noise_level) {
    wave noisy_wave = input_wave;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, noise_level);
    for (auto& sample : noisy_wave.y) {
        sample += d(gen);
    }
    matplot::ylabel("amplituda");
    matplot::title("Sygna³");
    matplot::plot(noisy_wave.x, noisy_wave.y)->color({ 1.0f, 0.08f, 0.58f });
    matplot::show();
    return noisy_wave;
}

void add_noise_to_audio(const std::string& inputAudioPath, const std::string& outputAudioPath, double noise_level) {
    AudioFile<double> audioFile;
    if (!audioFile.load(inputAudioPath)) {
        std::cerr << "Error loading audio file: " << inputAudioPath << std::endl;
        return;
    }

    int numSamples = audioFile.getNumSamplesPerChannel();
    int numChannels = audioFile.getNumChannels();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, noise_level);

    for (int channel = 0; channel < numChannels; ++channel) {
        for (int i = 0; i < numSamples; ++i) {
            audioFile.samples[channel][i] += d(gen);
        }
    }

    if (!audioFile.save(outputAudioPath)) {
        std::cerr << "Error saving audio file: " << outputAudioPath << std::endl;
    }
}

namespace py = pybind11;

PYBIND11_MODULE(cmake_example, m) {
    py::class_<wave>(m, "wave")
        .def(py::init<>())
        .def_readwrite("x", &wave::x)
        .def_readwrite("y", &wave::y);

    m.def("plot_sin", &plot_sinus, "sinus(x)");
    m.def("plot_cos", &plot_cosinus, "cosinus(x)");
    m.def("visualize", &visualize_audio, "visualisation");
    m.def("plot_sawtooth", &plot_sawtooth, "plot_sawtooth");
    m.def("plot_square", &plot_square, "square");
    m.def("sawtooth", &sawtooth, "sawtooth");
    m.def("square", &square, "square");
    m.def("sin", &sinus, "sinus(x)");
    m.def("cos", &cosinus, "cosinus(x)");
    m.def("dft", &plot_dft, "dft");
    m.def("idft", &plot_idft, "idft");
    m.def("DFT", &DFT, "DFT");
    m.def("noise", &add_noise, "noise");
    m.def("noise_audio", &add_noise_to_audio, "adding noise to audio");
}
