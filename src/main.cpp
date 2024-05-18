#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include "AudioFile.h"
#include <vector>
#include <cmath>

int add(int i, int j) {
    return i + j;
}

void sinus(double frequency) {

    std::vector<double> x;
    std::vector<double> y;

    for (int i = 0; i < 628; i++)
    {
        x.push_back(static_cast<double>(i) / 314);
        y.push_back(std::sin(3.14 * frequency * x[i]));
    }
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ sinusoidalny");
    matplot::plot(x, y)->color({ 1.0f, 0.08f, 0.58f });
    matplot::show();
}

void cosinus(double frequency) {
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i < 628; i++)
    {
        x.push_back(static_cast<double>(i) / 314);
        y.push_back(std::cos(3.14 * frequency * x[i]));
    }
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ cosinusoidalny");
    matplot::plot(x, y)->color({ 1.0f, 0.08f, 0.58f });;
    matplot::show();
}

void plot_sawtooth(int f, int fs, int num_periods) {
    double step = (2.0 * 32768) / (fs * f);
    std::vector<double> x;
    std::vector<double> sawtooth;
    double amplitude = 0;
    for (int p = 0; p < num_periods; ++p) {
        for (int i = 0; i < fs; ++i) {
            x.push_back(static_cast<double>(i + p * fs) / fs);
            sawtooth.push_back(amplitude);
            amplitude += step;
            if (amplitude >= 32768)
                amplitude -=65536 ;
        }
    }
    matplot::plot(x, sawtooth)->color({ 1.0f, 0.08f, 0.58f });
    matplot::xlabel("nr próbki");
    matplot::ylabel("amplituda");
    matplot::title("Sygna³ pi³okszta³tny");
    matplot::show();
}

void square(double amplitude, double period) {
    std::vector<double> x;
    std::vector<double> sawtooth;
    for (int i = 0; i < 628; i++) {
        x.push_back(static_cast<double>(i) / 314);
        sawtooth.push_back(amplitude * (2 * (x[i] / period - std::floor(x[i] / period + 0.5))));
    }
    std::vector<double> square;
    for (int i = 0; i < sawtooth.size(); i++) {
        square.push_back((sawtooth[i] >= 0) ? 1 : -1);
    }
    matplot::plot(x, square)->color({ 1.0f, 0.08f, 0.58f });
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

namespace py = pybind11;

PYBIND11_MODULE(cmake_example, m) {
    

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");
    m.def("sin", &sinus, "sinus(x)");
    m.def("cos", &cosinus, "cosinus(x)");
    m.def("visualize", &visualize_audio, "visualisation");
    m.def("plot_sawtooth", &plot_sawtooth, "plot_sawtooth");
    m.def("square", &square, "square");
}
