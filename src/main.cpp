#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
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
    std::vector<std::complex<double>> X;
    std::string name;
    wave(std::string name) : name(name), X(0) {}
    wave() : name("wave"), X(0){}
};

wave sinus(double frequency, double amplitude, double sampling_rate, double duration) {
    wave sinus("sin");
    int n_samples = static_cast<int>(sampling_rate * duration);
    for (int i = 0; i < n_samples; ++i) {
        double t = i / sampling_rate;
        sinus.x.push_back(t);
        sinus.y.push_back(amplitude * sin(2 * M_PI * frequency * t));
    }
    return sinus;
}
wave cosinus(double frequency, double amplitude, double sampling_rate, double duration) {
    wave cosine("cos");
    int n_samples = static_cast<int>(sampling_rate * duration);
    for (int i = 0; i < n_samples; ++i) {
        double t = i / sampling_rate;
        cosine.x.push_back(t);
        cosine.y.push_back(amplitude * cos(2 * M_PI * frequency * t));
    }
    return cosine;
}

wave sawtooth(double frequency, double amplitude, double sampling_rate, double duration) {
    wave sawtooth("sawtooth");
    int n_samples = static_cast<int>(sampling_rate * duration);
    for (int i = 0; i < n_samples; ++i) {
        double t = i / sampling_rate;
        sawtooth.x.push_back(t);
        sawtooth.y.push_back(amplitude * (2 * (t * frequency - floor(t * frequency + 0.5))));
    }
    return sawtooth;
}

wave square(double frequency, double amplitude, double sampling_rate, double duration) {
    wave square("square");
    int n_samples = static_cast<int>(sampling_rate * duration);
    for (int i = 0; i < n_samples; ++i) {
        double t = i / sampling_rate;
        square.x.push_back(t);
        square.y.push_back(amplitude * (sin(2 * M_PI * frequency * t) >= 0 ? 1 : -1));
    }
    return square;
}

wave DFT( wave input_wave) {
    wave output_wave = input_wave;
    output_wave.name += " after IDFT";
    int N = input_wave.y.size();
    output_wave.X.resize(N);
    output_wave.y.resize(N);
    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * M_PI * k * n / N;
            sum += input_wave.y[n] * std::complex<double>(std::cos(angle), std::sin(angle));
        }
        output_wave.X[k] = sum;
        output_wave.y[k] = std::abs(sum);
    }
    return output_wave;
}

wave IDFT( wave input_wave) {
    wave output_wave = input_wave;
    output_wave.name += " after IDFT";
    int N = input_wave.X.size();
    output_wave.y.resize(N);
    for (int n = 0; n < N; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k) {
            double angle = 2.0 * M_PI * k * n / N;
            sum += input_wave.X[k] * std::complex<double>(std::cos(angle), std::sin(angle));
        }
        output_wave.y[n] = std::real(sum) / N;
    }
    return output_wave;
}

void plot(wave input_wave) {
    matplot::title(input_wave.name);
    matplot::plot(input_wave.x, input_wave.y)->color({ 1.0f, 0.08f, 0.58f });
    matplot::show();
}


wave audio_to_wave( std::string audioFilePath) {
    AudioFile<double> audioFile;
    wave audio("audio");
    bool load = audioFile.load(audioFilePath);
    if (!load) {
        std::cerr << "Audio is not loaded" << std::endl;
        return audio;
    }
    int samples = audioFile.getNumSamplesPerChannel();
    int channels = audioFile.getNumChannels();
    if (samples <= 0) {
        std::cerr << "audio has no samples" << std::endl;
        return audio;
    }
    std::vector<double> audioData(audioFile.samples[0].begin(), audioFile.samples[0].end());
    std::vector<double> x(samples);
    for (int i = 0; i < samples; ++i) {
        x[i] = static_cast<double>(i) / audioFile.getSampleRate();
    }
    audio.x = x;
    audio.y = audioData;
    return audio;
}

wave add_noise(wave input_wave, double noise_level) {
    wave noisy_wave = input_wave;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, noise_level);
    for (auto& sample : noisy_wave.y) {
        sample += d(gen);
    }
    return noisy_wave;
}

void wave_to_audio(wave input_wave, std::string output_audio_path) {
    AudioFile<double> audio_file;
    audio_file.setSampleRate(44100);
    audio_file.setAudioBufferSize(1, input_wave.x.size());
    audio_file.samples[0] = input_wave.y;
    audio_file.save(output_audio_path);
}

namespace py = pybind11;

PYBIND11_MODULE(cmake_example, m) {
    py::class_<wave>(m, "wave")
        .def(py::init<>())
        .def_readwrite("x", &wave::x)
        .def_readwrite("y", &wave::y)
        .def_readwrite("name", &wave::name);

    m.def("audio_to_wave", &audio_to_wave, "audio to wave");
    m.def("sawtooth", &sawtooth, "sawtooth");
    m.def("square", &square, "square");
    m.def("sin", &sinus, "sinus(x)");
    m.def("cos", &cosinus, "cosinus(x)");
    m.def("DFT", &DFT, "DFT");
    m.def("IDFT", &IDFT, "IDFT");
    m.def("noise", &add_noise, "noise");
    m.def("wave_to_audio", &wave_to_audio, "wave to audio");
    m.def("plot", &plot, "plot");
}
