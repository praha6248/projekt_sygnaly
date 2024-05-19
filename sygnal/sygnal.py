import sys
sys.path.append('../out/build/x64-debug')
import cmake_example as lib
#lib.wave_to_audio(lib.noise(lib.visualize("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav"), 0.02), "C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/noise.wav");
#lib.noise_audio("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav", "C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/noise.wav", 0.02)
lib.plot(lib.square(1, 45, 5))
lib.plot(lib.DFT(lib.square(1, 45, 5)))
lib.plot(lib.IDFT(lib.DFT(lib.square(1, 45, 5))))
