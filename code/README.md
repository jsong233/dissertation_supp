## Compile

### compile on `logopolis` and `darth`:

`clang -O3 -stdlib=libc++ -std=c++17 -lc++ cellclass_def.cpp main-file-name.cpp -I/Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers -lsundials_cvode -lsundials_nvecserial -lblas -llapack -lm -o executable-file-name`

### compile on MacOS (M1):

`arch -x86_64 clang -O3 -stdlib=libc++ -std=c++17 -lc++ cellclass_def.cpp main-file-name.cpp -I/Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers -lsundials_cvode -lsundials_nvecserial -lblas -llapack -lm -o executable-file-name`

### compile on `castrovalva` and `sontar`

`clang -O3 -stdlib=libc++ -std=c++11 -lc++ cellclass_def.cpp main-file-name.cpp -I/usr/local/include -L/usr/local/lib -I/Library/Developer/CommandLineTools/SDKs/MacOSX10.14.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers -lsundials_cvode -lsundials_nvecserial -lblas -llapack -lm -o Mar16`

## Execute

`nohup ./executable-file-name > temp.dat &`
