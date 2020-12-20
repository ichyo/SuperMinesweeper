import subprocess

def run(d, r, case_num=1000, n=50):
    m = int(n * n * r / 100)
    print('n={}'.format(n), flush=True)
    print('m={}'.format(m), flush=True)
    print('d={}'.format(d), flush=True)
    result = subprocess.run(['./run_test.sh', '-l', str(n), '-d', str(d), '-m', str(m), '-n', str(case_num)], check=True, stdout=subprocess.PIPE, encoding='utf-8').stdout
    result = str(result)
    csv_path = result.splitlines()[-1].split(" ")[-1].strip()
    print(csv_path)
    subprocess.run(['python', './mine_hit_stats.py', csv_path], check=True)

ds = [1, 2, 4, 5, 8, 9, 10]
for d in [1]:
    for r in range(10, 31, 2):
        run(d, r)
