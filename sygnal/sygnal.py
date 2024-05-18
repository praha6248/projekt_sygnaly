import sys
sys.path.append('../out/build/x64-debug')
import cmake_example as lib
#lib.visualize("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav")
#lib.plot_sin(5)
#lib.plot_square(1, 45, 5)
sine_wave=lib.square(1,45,5)
lib.dft(sine_wave)
sine_wave=lib.DFT(sine_wave)
lib.idft(sine_wave)