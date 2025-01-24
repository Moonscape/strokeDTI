import os
import urllib.request
import tarfile

# URL of the trained model archive
MODEL_URL = "https://github.com/Moonscape/strokeDTI/releases/download/v0.1.0-model/trained_model.tar.gz"


def main():
    # Define the path to the 'data' directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, "data")

    # Define the path to the 'trained_model' directory
    model_dir = os.path.join(data_dir, "trained_model")

    # Check if the 'trained_model' directory already exists
    if os.path.isdir(model_dir):
        print("Model already downloaded.")
        return
    else:
        print("Downloading and extracting the model...")

    # Ensure the 'data' directory exists
    os.makedirs(data_dir, exist_ok=True)

    # Path where the model archive will be saved
    model_archive_path = os.path.join(data_dir, "trained_model.tar.gz")

    try:
        # Download the model archive
        print(f"Downloading model from {MODEL_URL}...")
        urllib.request.urlretrieve(MODEL_URL, model_archive_path)
        print("Download completed.")

        # Extract the archive
        print("Extracting the model...")
        with tarfile.open(model_archive_path, "r:gz") as tar:
            tar.extractall(path=data_dir)
        print("Extraction completed.")

    except Exception as e:
        print(f"An error occurred: {e}")
        return

    finally:
        # Optionally, remove the archive file
        if os.path.exists(model_archive_path):
            os.remove(model_archive_path)
            print("Cleaned up the archive file.")

    print("Model is ready to use.")


if __name__ == "__main__":
    main()
