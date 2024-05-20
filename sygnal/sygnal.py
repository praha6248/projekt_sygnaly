import sys
sys.path.append('../out/build/x64-debug')
import cmake_example as lib
#lib.wave_to_audio(lib.noise(lib.audio_to_wave("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav"), 0.1), "C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/noise.wav")
#lib.noise_audio("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav", "C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/noise.wav", 0.02)
wave =lib.sin(3,2,100,2)
lib.plot(wave)
dft_wave=lib.DFT(wave)
lib.plot(dft_wave)
lib.plot(lib.IDFT(wave))