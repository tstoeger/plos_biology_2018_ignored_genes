import glob
import os.path
import argparse
import subprocess
import time
import numpy as np


def submit_zgbrh_stack_all(datasets, p, e, n, a, r):

    maximal_jobs_to_submit = 450

    minimally_required_predictions = r
    base_dir = datasets
    print('the basedir is', base_dir)

    # consider making testing for folders more explicit
    projects = glob.glob(os.path.join(base_dir, '1*'))
    print('------found projects ------')
    for dummy in sorted(projects):
        print(dummy)

    print('------------')

    currently_sumbitted = 0
    for i in projects:
        # print('entering project', i)

        _, dataset = os.path.split(i)

        i_scripts = os.path.join(i, 'scripts')

        estimators = int(np.round(e))
        n = int(n)
        predictor_name = 'zgbrh_p{}_e{}'.format(p, estimators)

        submit_jobs = False
        p_folder = os.path.join(i, predictor_name)

        # print('anticipating results folder', p_folder)

        existing = glob.glob(os.path.join(
            p_folder, 'prediction_*'))
        amount_of_files = len(existing)

        print('project {} has {} pre-existing results'.format(
            i, amount_of_files))

        if amount_of_files < minimally_required_predictions:
            # print(amount_of_files, 'is smaller than',
            #       minimally_required_predictions)
            submit_jobs = True
        else:
            submit_jobs = False
            print(
                'Will not submit, as there are already {} files, and only {} are requested.'.format(
                    amount_of_files, minimally_required_predictions))

        for job in np.arange(a):
            if currently_sumbitted >= maximal_jobs_to_submit:
                submit_jobs = False
                print(
                    'Can not submit. Already submitted {} jobs.'.format(
                        maximal_jobs_to_submit))

            if submit_jobs:
                print('Going to submit {} as job number {}'.format(
                    p_folder, currently_sumbitted))
                currently_sumbitted = currently_sumbitted + 1

                shellscript = "{}/{}_{}.sh".format(
                    i_scripts,
                    predictor_name,
                    str(int(round(time.time() * 1000)))
                )

                subprocess.call(["mkdir", "-p", i_scripts])
                with open(shellscript, 'w') as f:
                    f.write("#!/bin/bash \n")
                    f.write("#MSUB -q normal \n")
                    f.write("#MSUB -l walltime=48:00:00 \n")
                    f.write("#MSUB -M thomas.stoeger@northwestern.edu \n")
                    f.write("#MSUB -j oe \n")
                    f.write("#MSUB -N %s/ \n" % dataset)
                    f.write("#MSUB -l nodes=1:ppn=4 \n")
                    f.write("module load python/anaconda3.6 \n")
                    f.write("cd {} \n".format(base_dir))
                    f.write("source activate fame_predictions \n")
                    f.write(
                        "python ./predict_via_zgbrh_importance.py -i {} -p {} -e {} -n {} \n".format(
                            i, p, e, n))
                subprocess.call(["msub", "%s" % shellscript])

                time.sleep(10)   # so not to overwhelm quest's scheduler


def main():
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Submit zgbrh for a given prediction set.')
    parser.add_argument('-d', type=str,
                        help='name/of/datasets')
    parser.add_argument('-n', type=int,
                        help='number/of/iterations/per/job')
    parser.add_argument('-p', type=int,
                        help='percentiles/to/train')
    parser.add_argument('-e', type=int,
                        help='estimators/for/model')
    parser.add_argument('-a', type=int,
                        help='amount/of/jobs')
    parser.add_argument('-r', type=int,
                        help='required/predictions')

    args = parser.parse_args()

    submit_zgbrh_stack_all(
        args.d,    # datasets
        args.p,    # percentiles to train
        args.e,    # estimators
        args.n,    # number of iterations
        args.a,    # amount of jobs to be submitted
        args.r,    # amount of required predictions
    )
