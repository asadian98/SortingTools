
import sys
import deeplabcut

def analyze(path_config_file, videofile_path):

    deeplabcut.analyze_videos(path_config_file,videofile_path, videotype='.mp4',save_as_csv=True)

    deeplabcut.filterpredictions(path_config_file, videofile_path, save_as_csv=True)


if __name__ == "__main__":
    path_config_file = sys.argv[1]  # First argument from the command line
    videofile_path = sys.argv[2:]   # Remaining arguments as a list for video files
    analyze(path_config_file, videofile_path)