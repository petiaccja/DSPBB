add_executable(01_Basics "01_Basics.cpp")
add_executable(02_FFT "02_FFT.cpp")
add_executable(03_IIR_Filtering "03_IIR_Filtering.cpp" "RtAudio/RtAudio.cpp" "Tools/PlaySound.cpp")
add_executable(04_FIR_Filtering "04_FIR_Filtering.cpp" "RtAudio/RtAudio.cpp" "Tools/LoadSound.cpp" "Tools/PlaySound.cpp")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR}/examples)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/examples)

conan_cmake_configure(
	REQUIRES
		libsndfile/1.0.31
	GENERATORS cmake_find_package
)
conan_cmake_install(PATH_OR_REFERENCE . BUILD missing SETTINGS ${conan_settings} OUTPUT_QUIET)

if (${CMAKE_CXX_COMPILER} MATCHES "MSVC")
	target_compile_options(01_Basics PRIVATE "/bigobj")
	target_compile_options(02_FFT PRIVATE "/bigobj")
	target_compile_options(03_IIR_Filtering PRIVATE "/bigobj")
	target_compile_options(04_FIR_Filtering PRIVATE "/bigobj")
endif()

find_package(SndFile REQUIRED)

target_link_libraries(01_Basics DSPBB)
target_link_libraries(02_FFT DSPBB)
target_link_libraries(03_IIR_Filtering DSPBB)
target_link_libraries(04_FIR_Filtering DSPBB SndFile::SndFile)

set(DSPBB_EXAMPLES_DATA "\"${CMAKE_CURRENT_SOURCE_DIR}/data/\"")
target_compile_definitions(02_FFT PRIVATE "-DDSPBB_EXAMPLES_DATA=${DSPBB_EXAMPLES_DATA}")
target_compile_definitions(04_FIR_Filtering PRIVATE "-DDSPBB_EXAMPLES_DATA=${DSPBB_EXAMPLES_DATA}")

if (${WIN32})
	set(RTAUDIO_API "__WINDOWS_WASAPI__")
elseif (${APPLE})
	set(RTAUDIO_API "__MACOSX_CORE__")
elseif (${UNIX})
	find_package(ALSA REQUIRED)
	if (${ALSA_FOUND})
		set(RTAUDIO_API "__LINUX_ALSA__")
		target_link_libraries(03_IIR_Filtering ALSA::ALSA)
		target_link_libraries(04_FIR_Filtering ALSA::ALSA)
	else()
		set(RTAUDIO_API "__RTAUDIO_DUMMY__")
	endif()
else()
	set(RTAUDIO_API "__RTAUDIO_DUMMY__")
endif()

target_compile_definitions(03_IIR_Filtering PRIVATE "-D${RTAUDIO_API}=1")
target_compile_definitions(04_FIR_Filtering PRIVATE "-D${RTAUDIO_API}=1")