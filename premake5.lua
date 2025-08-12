workspace "Master-Equation"
configurations {"Debug", "Release"}
files{"vendor/Lua/src/*.c"}
removefiles{"vendor/Lua/src/lua.c", "vendor/Lua/src/luac.c"}

language "C++"
cppdialect "C++20"
toolset "clang"
buildoptions {"-Wall", "-march=native"}
links {}
-- "-fopenmp for compiler", "omp" for linker
targetdir "bin"

project "main"
kind "ConsoleApp"
files {"main/src/*.cpp", "tools/*.cpp", "math/*.cpp"}
removefiles {}
libdirs {}
links {}
-- linkoptions{"-framework Accelerate", "/opt/local/lib/lapack/liblapacke.dylib"} -- Apple specific
location "main/build"
-- buildoptions{"-DACCELERATE_NEW_LAPACK"} -- Apple specific

filter {"configurations:Debug"}
	buildoptions {"-DDEBUG", "-g"}
	targetname "meq-d"

filter {"configurations:Release"}
	buildoptions {"-DRELEASE", "-O3", "-Wno-deprecated-anon-enum-enum-conversion"}
	targetname "meq"


project "occupations"
kind "ConsoleApp"
files {"occupations/src/*.cpp", "tools/*.cpp", "math/*.cpp"}
removefiles {}
libdirs {}
links {}
-- linkoptions{"-framework Accelerate", "/opt/local/lib/lapack/liblapacke.dylib"} -- Apple specific
location "main/build"
-- buildoptions{"-DACCELERATE_NEW_LAPACK"} -- Apple specific

filter {"configurations:Debug"}
	buildoptions {"-DDEBUG", "-g"}
	targetname "occ-d"

filter {"configurations:Release"}
	buildoptions {"-DRELEASE", "-O3", "-Wno-deprecated-anon-enum-enum-conversion"}
	targetname "occ"


