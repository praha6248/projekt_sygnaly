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
    std::string name;
    wave(const std::string& name) : name(name) {}
    wave() : name("wave") {}
};

wave sinus(double frequency) {
    wave sinus("sin");
    for (int i = 0; i < 628; i++)
    {
        sinus.x.push_back(static_cast<double>(i) / 314);
        sinus.y.push_back(std::sin(3.14 * frequency * sinus.x[i]));
    }
    return sinus;
}
wave cosinus(double frequency) {
    wave cosinus("cos");
    for (int i = 0; i < 628; i++)
    {
        cosinus.x.push_back(static_cast<double>(i) / 314);
        cosinus.y.push_back(std::cos(3.14 * frequency * cosinus.x[i]));
    }
    return cosinus;
}

wave sawtooth(int f, int fs, int num_periods) {
    wave sawtooth("sawtooth");
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
    wave square("square");
    for (int p = 0; p < num_periods; ++p) {
        for (int i = 0; i < fs; ++i) {
            square.x.push_back(static_cast<double>(i + p * fs) / fs);
            square.y.push_back(i < fs / 2 ? -32768 : 32768);
        }
    }
    return square;
}

wave DFT(const wave& input_wave) {
    int N = input_wave.y.size();
    wave output_wave("DFT"+input_wave.name);
    output_wave.x = input_wave.x;
    output_wave.y.resize(N);
    std::vector<std::complex<double>> X(N);
    for (int k = 0; k < N; ++k) {
        std::complex<double> sum = 0.0;
        for (int n = 0; n < N; ++n) {
            double angle = -2 * M_PI * k * n / N;
            sum += input_wave.y[n] * std::exp(std::complex<double>(0, angle));
        }
        X[k] = sum;
        output_wave.y[k] = std::abs(sum);
    }
    return output_wave;
}

wave IDFT(const wave& input_wave) {
    int N = input_wave.y.size();
    wave output_wave("IDFT");
    output_wave.x.resize(N);
    output_wave.y.resize(N);

    for (int n = 0; n < N; ++n) {
        std::complex<double> sum = 0.0;
        for (int k = 0; k < N; ++k) {
            double angle = 2 * M_PI * k * n / N;
            sum += std::complex<double>(input_wave.y[k]) * std::exp(std::complex<double>(0, angle));
        }
        output_wave.x[n] = input_wave.x[n];
        output_wave.y[n] = std::real(sum) / N;
    }

    return output_wave;
}

void plot(wave input_wave) {
    matplot::ylabel("amplituda");
    matplot::title(input_wave.name);
    matplot::plot(input_wave.x, input_wave.y)->color({ 1.0f, 0.08f, 0.58f });
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
        .def_readwrite("y", &wave::y)
        .def_readwrite("name", &wave::name);

    m.def("visualize", &visualize_audio, "visualisation");
    m.def("sawtooth", &sawtooth, "sawtooth");
    m.def("square", &square, "square");
    m.def("sin", &sinus, "sinus(x)");
    m.def("cos", &cosinus, "cosinus(x)");
    m.def("DFT", &DFT, "DFT");
    m.def("IDFT", &IDFT, "IDFT");
    m.def("noise", &add_noise, "noise");
    m.def("noise_audio", &add_noise_to_audio, "adding noise to audio");
}
