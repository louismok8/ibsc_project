import sys
from pathlib import Path

PHASE_ROOT = Path(__file__).resolve().parent
if str(PHASE_ROOT) not in sys.path:
    sys.path.insert(0, str(PHASE_ROOT))

from step1_load_data import load_full_dataset


def main():
    dataset = load_full_dataset()
    dataset.verify_all_geometry()
    print("✓ Step 2 complete: all geometry valid.")


if __name__ == "__main__":
    main()
