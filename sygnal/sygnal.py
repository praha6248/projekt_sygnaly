import sys
sys.path.append('../out/build/x64-debug')
import cmake_example as lib
#lib.visualize("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav")
#lib.plot_sin(5)
#lib.plot_square(1, 45, 5)
#sine_wave=lib.sin(3)
#lib.plot_sin(3)
#lib.dft(sine_wave)
#sine_wave=lib.DFT(sine_wave)
#lib.idft(sine_wave)
#lib.noise(sine_wave, 0.5)
lib.noise_audio("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav", "C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/noise.wav", 0.02)
