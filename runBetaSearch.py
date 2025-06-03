#!/usr/bin/env python3

import subprocess
import os
from pathlib import Path
import time

def run_deeptmhmm(file_path, max_retries=3):
    base = file_path.stem
    result_dirname = f"{base}_biolib"

    # Skip if already processed
    if Path(result_dirname).is_dir():
        print(f"Skipping {file_path.name}: results already exist in {result_dirname}")
        return

    for attempt in range(1, max_retries + 1):
        print(f"Running DeepTMHMM on {file_path.name}, attempt {attempt}...")

        try:
            result = subprocess.run(
                ["biolib", "run", "--local", "DTU/DeepTMHMM:1.0.24", "--fasta", str(file_path)],
                stdout=None,              # show live stdout
                stderr=subprocess.PIPE,   # capture stderr for checking errors
                text=True
            )

            # Print stderr in all cases
            if result.stderr:
                print(result.stderr)

            # Retry only on traceback errors
            if "Traceback (most recent call last):" in result.stderr:
                print(f"Server-side error detected on attempt {attempt} for {base}. Retrying...")
                time.sleep(5)
                continue

        except Exception as e:
            print(f"Exception occurred running biolib on {base}: {e}")
            time.sleep(5)
            continue

        # Rename result folder if it exists
        temp_result_dir = Path("biolib_results")
        if temp_result_dir.is_dir():
            os.rename(temp_result_dir, result_dirname)
            print(f"Renamed biolib_results → {result_dirname}")
        else:
            print(f"Completed run for {base}, but no biolib_results folder found.")
        break

def main():
    for file in Path(".").glob("*.faa"):
        run_deeptmhmm(file)

if __name__ == "__main__":
    main()
