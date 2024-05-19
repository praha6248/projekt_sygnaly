import sys
sys.path.append('../out/build/x64-debug')
import cmake_example as lib
lib.visualize("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav")
lib.noise_audio("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav", "C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/noise.wav", 0.02)
