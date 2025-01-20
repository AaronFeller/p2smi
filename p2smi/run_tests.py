import subprocess
import sys

def run_tests():
    test_files = [
        "test_conversion.py",
        "test_modify_smiles.py",
        "test_smilesgen.py",
        "test_synthrules.py",  # Add all your test files here
    ]

    for test_file in test_files:
        print(f"Running {test_file}...")
        result = subprocess.run([sys.executable, test_file], capture_output=True, text=True)

        if result.returncode == 0:
            print(f"{test_file} passed.\n")
        else:
            print(f"{test_file} failed.\n")
            print("Output:\n", result.stdout)
            print("Errors:\n", result.stderr)

    print("All tests completed.")

if __name__ == "__main__":
    run_tests()