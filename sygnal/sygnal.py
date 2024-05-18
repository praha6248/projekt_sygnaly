import sys
sys.path.append('../out/build/x64-debug')
import cmake_example as lib
print(lib.add(2, 4))
#lib.visualize("C:/Users/natal/OneDrive/Pulpit/projektsygnaly/cmake_example/test-audio.wav")
lib.plot_sawtooth(1, 48, 5)
lib.cos(5)
#lib.square(10, 5)