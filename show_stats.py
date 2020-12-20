import argparse
import os.path
import yaml
import csv
import collections
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

meta = {}
for seed in seeds:
    input_path = os.path.join(base_dir, 'io', '{}.err'.format(seed))
    reason = ""
    runtime = 0
    uncover_ratio = 0.0
    mine_hit = 0
    guess_count = 0
    random_guess_count = 0
    with open(input_path) as f:
        for s in f.readlines():
            if s.startswith("reason: "):
                reason = s.split(" ")[1].strip()
            if s.startswith("runtime: "):
                runtime = int(s.split(" ")[1].strip())
            if s.startswith("uncover_ratio: "):
                uncover_ratio = float(s.split(" ")[1].strip())
            if s.startswith("mine_hit: "):
                mine_hit = int(s.split(" ")[1].strip())
            if s.startswith("guess_count: "):
                guess_count = int(s.split(" ")[1].strip())
            if s.startswith("random_guess_count: "):
                random_guess_count = int(s.split(" ")[1].strip())
        if reason == "":
            print('Warning: seed={} has no reason'.format(seed))
    meta[seed] = {
        'reason': reason,
        'runtime': runtime,
        'uncover_ratio': uncover_ratio,
        'mine_hit': mine_hit,
        'guess_count': guess_count,
        'random_guess_count': random_guess_count,
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

print('=== param stats ===')
ds = sorted(set(p['d'] for p in params.values()))
for d in ds:
    print('d={}'.format(d))
    values = score_values[[seed - first_seed for seed in seeds if params[seed]['d'] == d]]
    print('  cnt: {}'.format(len(values)))
    print('  avg: {:.4f}'.format(values.mean()))
    print('  25t: {:.4f}'.format(np.percentile(values, 25)))
    print('  50t: {:.4f}'.format(np.median(values)))
    print('  75t: {:.4f}'.format(np.percentile(values, 75)))

print('=== meta stats ===')
reasons = sorted(set(m['reason'] for m in meta.values()))
for reason in reasons:
    print('reason={}'.format(reason))
    values = score_values[[seed - first_seed for seed in seeds if meta[seed]['reason'] == reason]]
    print('  cnt: {}'.format(len(values)))
    print('  avg: {:.4f}'.format(values.mean()))
    print('  25t: {:.4f}'.format(np.percentile(values, 25)))
    print('  50t: {:.4f}'.format(np.median(values)))
    print('  75t: {:.4f}'.format(np.percentile(values, 75)))



csv_path = os.path.join(base_dir, 'scores.csv')
with open(csv_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['seed', 'n', 'm', 'd', 'density', 'score', 'reason', 'runtime', 'uncover_ratio', 'mine_hit', 'guess_count', 'random_guess_count'])
    for seed in seeds:
        n = params[seed]['n']
        m = params[seed]['m']
        d = params[seed]['d']
        reason = meta[seed]['reason']
        runtime = meta[seed]['runtime']
        uncover_ratio = meta[seed]['uncover_ratio']
        mine_hit = meta[seed]['mine_hit']
        guess_count = meta[seed]['guess_count']
        random_guess_count = meta[seed]['random_guess_count']
        writer.writerow([seed, n, m, d, float(m) / float(n * n), scores[seed], reason, runtime, uncover_ratio, mine_hit, guess_count, random_guess_count])
print('Written into {}'.format(csv_path))
