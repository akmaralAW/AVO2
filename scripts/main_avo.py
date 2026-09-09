import subprocess
from pathlib import Path
import sys
import os
import yaml
sys.path.append(os.getcwd())

def run_avo(filename_env, folder, timelimit):
    
    filename_result = Path(folder) / "result_avo.yaml"
    filename_stats = Path(folder) / "stats.yaml"

    cmd = ["./planners/AVO2/build/AVOPlanner", 
        filename_env,
        filename_result,
        filename_stats,
        ]
    
    print(subprocess.list2cmdline(cmd))
    try:
        with open("{}/log.txt".format(folder), 'w') as logfile:
          result = subprocess.run(cmd, timeout=timelimit, stdout=logfile, stderr=logfile)

    except subprocess.TimeoutExpired as e:
        print(f"Command timed out after {timelimit} seconds.")
        cmd_str = " ".join(str(arg) for arg in cmd)
        subprocess.run(["pkill", "-f", " ".join(cmd_str)])
        # Check if stats file exists and load stats
        filename_stats = f"{folder}/stats.yaml"
        try:
            with open(filename_stats, 'r') as file:
                stats = yaml.safe_load(file)
                if stats and "stats" in stats and len(stats["stats"]):
                    print("AVO succeeded!")
                else:
                    print("AVO failed!")
        except FileNotFoundError:
            print(f"Stats file {filename_stats} not found.")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

