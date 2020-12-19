import argparse
import os.path
import yaml
import csv
import numpy as np

parser = argparse.ArgumentParser(description='Show stats of test')

parser.add_argument('base_dir')

args = parser.parse_args()

base_dir = args.base_dir
if not os.path.isdir(base_dir):
    raise Exception('{} is not directory'.format(base_dir))

info_file = os.path.join(base_dir, 'info.yml')
if not os.path.isfile(info_file):
    raise Exception('{} is not found'.format(info_file))

with open(info_file) as f:
    info = yaml.safe_load(f)

first_seed = info['first_seed']
case_num = info['case_num']
seeds = list(range(first_seed, first_seed + case_num))

scores_path = os.path.join(base_dir, 'scores.txt')
scores = {}
with open(scores_path) as f:
    for line in f.readlines():
        seed, score = line.split('=')
        seed = int(seed)
        score = float(score)
        if score < 0.0:
            score = 0.0 # convert -1 to 0 for stats
        scores[seed] = score

params = {}
for seed in seeds:
    input_path = os.path.join(base_dir, 'io', '{}.in'.format(seed))
    with open(input_path) as f:
        param_n = int(f.readline())
        param_m = int(f.readline())
        param_d = int(f.readline())
    params[seed] = {
        'n': param_n,
        'm': param_m,
        'd': param_d,
    }

score_values = np.array(list(scores.values()))

print('=== score stats ===')
print('avg: {:.4f}'.format(score_values.mean()))
print('std: {:.4f}'.format(score_values.std()))
print('min: {:.4f}'.format(score_values.min()))
print('05t: {:.4f}'.format(np.percentile(score_values, 5)))
print('25t: {:.4f}'.format(np.percentile(score_values, 25)))
print('50t: {:.4f}'.format(np.median(score_values)))
print('75t: {:.4f}'.format(np.percentile(score_values, 75)))
print('95t: {:.4f}'.format(np.percentile(score_values, 95)))

csv_path = os.path.join(base_dir, 'scores.csv')
with open(csv_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['seed', 'n', 'm', 'd', 'density', 'score'])
    for seed in seeds:
        n = params[seed]['n']
        m = params[seed]['m']
        d = params[seed]['d']
        writer.writerow([seed, n, m, d, float(m) / float(n * n), scores[seed]])
print('Written into {}'.format(csv_path))
