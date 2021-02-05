#define CATCH_CONFIG_RUNNER
#include <Catch2/catch.hpp>


int main(int argc, char* argv[]) {
	int result = Catch::Session().run(argc, argv);

	return result;
}
