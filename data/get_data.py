import os
import requests

instances = {
    "ta15x15": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai15_15.txt",
    "ta20x15": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai20_15.txt",
    "ta20x20": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai20_20.txt",
    "ta30x15": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai30_15.txt",
    "ta30x20": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai30_20.txt",
    "ta50x15": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai50_15.txt",
    "ta50x20": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai50_20.txt",
    "ta100x20": "http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/jobshop.dir/tai100_20.txt",
}


def parse_instance_file(text, out_dir):
    lines = text.strip().splitlines()
    i = 0
    instance_id = 1

    while i < len(lines):
        # Skip lines until we find a line with exactly 6 numbers (the instance header)
        while i < len(lines):
            tokens = lines[i].strip().split()
            if len(tokens) == 6 and all(t.isdigit() for t in tokens):
                num_jobs = int(tokens[0])
                num_machines = int(tokens[1])
                i += 1
                break
            i += 1
        else:
            break  # No more instances

        # Find "Times" section
        while i < len(lines) and not lines[i].strip().lower().startswith("times"):
            i += 1
        i += 1  # Skip the "Times" line

        processing_times = []
        for _ in range(num_jobs):
            processing_times.append(lines[i].strip())
            i += 1

        # Find "Machines" section
        while i < len(lines) and not lines[i].strip().lower().startswith("machines"):
            i += 1
        i += 1  # Skip the "Machines" line

        machines = []
        for _ in range(num_jobs):
            row = list(map(int, lines[i].strip().split()))
            zero_indexed = [x - 1 for x in row]
            machines.append(" ".join(map(str, zero_indexed)))
            i += 1

        # Save the extracted data
        with open(os.path.join(out_dir, f"processing_times_{instance_id}.txt"), "w") as f:
            f.write("\n".join(processing_times))

        with open(os.path.join(out_dir, f"machines_{instance_id}.txt"), "w") as f:
            f.write("\n".join(machines))

        instance_id += 1


def main():
    for name, url in instances.items():
        print(f"Processing {name}...")
        os.makedirs(name, exist_ok=True)
        response = requests.get(url)
        if response.status_code == 200:
            parse_instance_file(response.text, name)
        else:
            print(f"Failed to download {url}")

if __name__ == "__main__":
    main()