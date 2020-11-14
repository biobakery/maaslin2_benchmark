'''Run with $python model_comparison_workflow.py --threads (# of cores) --time (time in minutes per job)
--memory (in megabytes per task) --grid-jobs (# of jobs to run at a time --parameters (txt file containing parameters)
--input (input directory)  --output (output directory) [these two are the input and output for ANADAMA2, not the rscript]
(also change working directory in the data_generation_parameters.txt)
'''
from shared_workflow import workflow
from shared_workflow import param_list_evaluation
from shared_workflow import args
from shared_workflow import working_directory

for param in param_list_evaluation:
    if param[7]:
        output = 'output_filename'

        workflow.add_task_gridable('Rscript Script_name.R' + ' --metadataType ' + param[6] + ' --nMicrobes '
        + param[0] + ' --nSamples ' + param[1] + ' --spikePerc ' + param[2] + ' --spikeStrength ' + param[3] + ' --nMetadata '
        + param[4] + ' --spikeMeta ' + param[5]  + ' --noZeroInflate' + ' --methodName ' + param[-1]
        ' --workingDirectory' + ' ' + str(working_directory),time=args.time,mem=args.memory,cores=args.threads,target=output)

    else:
        output = 'output_filename'

        workflow.add_task_gridable('Rscript Script_name.R' + ' --metadataType ' + param[6] + ' --nMicrobes '
        + param[0] + ' --nSamples ' + param[1] + ' --spikePerc ' + param[2] + ' --spikeStrength ' + param[3] + ' --nMetadata '
        + param[4] + ' --spikeMeta ' + param[5] + ' --methodName ' + param[-1] + ' --workingDirectory' + ' ' + str(working_directory)
        ,time=args.time,mem=args.memory,cores=args.threads,target=output)


workflow.go()
