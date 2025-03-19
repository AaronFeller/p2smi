import subprocess
import sys


def run_tests():
    test_files = [
        "tests/test_chemMods.py",
        "tests/test_chemProps.py",
        "tests/test_fasta2smi.py",
        "tests/test_genPeps.py",
        "tests/test_synthRules.py",
        # Add more test files as needed
    ]

    for test_file in test_files:
        print(f"Running {test_file}...\n")
        result = subprocess.run(
            [sys.executable, "-m", "pytest", test_file, "-v"],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            print(f"{test_file} passed.\n")
        else:
            print(f"{test_file} failed.\n")
            print("Output:\n", result.stdout)
            print("Errors:\n", result.stderr)

    print("All tests completed.")


if __name__ == "__main__":
    run_tests()
