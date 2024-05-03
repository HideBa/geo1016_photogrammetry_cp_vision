#include "calibration.h"
#include <easy3d/fileio/resources.h>
#include <easy3d/util/logging.h>


using namespace easy3d;

int main(int argc, char** argv) {
    // the model file.
    const std::string model_file = resource::directory() + "/data/corner.obj";

    try {
        Calibration viewer("Calibration", model_file);

        // Run the viewer
        viewer.run();
    } catch (const std::runtime_error &e) {
        LOG(ERROR) << "caught a fatal error: " + std::string(e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
