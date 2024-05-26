#include "opencv2/opencv.hpp"

int main() {
    const std::string input = std::string(DATA_DIR) + "blais.mp4";
    const std::string cover = std::string(DATA_DIR) + "blais.jpg";
    size_t min_inlier_num = 100;

    // Load the object image and extract features
    cv::Mat obj_image = cv::imread(cover);
    if (obj_image.empty())
        return -1;

    cv::Ptr<cv::FeatureDetector> fdetector = cv::ORB::create();
    cv::Ptr<cv::DescriptorMatcher> fmatcher = cv::DescriptorMatcher::create("BruteForce-Hamming");
    std::vector<cv::KeyPoint> obj_keypoint;
    cv::Mat obj_descriptor;
    fdetector->detectAndCompute(obj_image, cv::Mat(), obj_keypoint, obj_descriptor);
    if (obj_keypoint.empty() || obj_descriptor.empty())
        return -1;
    fmatcher->add(obj_descriptor);

    // Open a video
    cv::VideoCapture video;
    if (!video.open(input))
        return -1;

    // Run pose estimation and camera calibration together
    while (true) {
        // Grab an image from the video
        cv::Mat image;
        video >> image;
        if (image.empty())
            break;

        // Extract features and match them to the object features
        std::vector<cv::KeyPoint> img_keypoint;
        cv::Mat img_descriptor;
        fdetector->detectAndCompute(image, cv::Mat(), img_keypoint, img_descriptor);
        if (img_keypoint.empty() || img_descriptor.empty())
            continue;
        std::vector<cv::DMatch> match;
        fmatcher->match(img_descriptor, match);
        if (match.size() < min_inlier_num)
            continue;
        std::vector<cv::Point3f> obj_points;
        std::vector<cv::Point2f> obj_project, img_points;
        for (auto m = match.begin(); m < match.end(); m++) {
            obj_points.push_back(cv::Point3f(obj_keypoint[m->trainIdx].pt));
            obj_project.push_back(obj_keypoint[m->trainIdx].pt);
            img_points.push_back(img_keypoint[m->queryIdx].pt);
        }

        // Determine whether each matched feature is an inlier or not
        cv::Mat inlier_mask = cv::Mat::zeros(int(match.size()), 1, CV_8U);
        cv::Mat H = cv::findHomography(img_points, obj_project, inlier_mask, cv::RANSAC, 2);
        cv::Mat image_result;
        cv::drawMatches(image, img_keypoint, obj_image, obj_keypoint, match, image_result, cv::Vec3b(0, 0, 255),
                        cv::Vec3b(0, 127, 0), inlier_mask);

        size_t inlier_num = static_cast<size_t>(cv::sum(inlier_mask)[0]);

        // Show the image
        cv::String info = cv::format("Inliers: %d (%d%%)", inlier_num, 100 * inlier_num / match.size());
        cv::putText(image_result, info, cv::Point(5, 15), cv::FONT_HERSHEY_PLAIN, 1, cv::Vec3b(0, 255, 0));
        cv::imshow("Real-time Image Matching", image_result);
        int key = cv::waitKey(1);
        if (key == 27)  // 'ESC' key: Exit
            break;
    }

    video.release();
    return 0;
}
