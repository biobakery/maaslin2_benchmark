'''Run with $python data_generation_workflow.py --threads (# of cores) --time (time in minutes per job)
--memory (in megabytes per task) --grid-jobs (# of jobs to run at a time --parameters (txt file containing parameters)
--input (input directory)  --output (output directory) [these two are the input and output for ANADAMA2, not the rscript]
(also change input and output directory in the data_generation_parameters.txt)
'''

from shared_workflow import workflow
from shared_workflow import param_list_generation
from shared_workflow import args
from shared_workflow import working_directory 

for param in param_list_generation:
    if param[-1]:
        output = 'noZeroInflate' + '_' + param[6] + '_' + param[0] + '_' + param[1] + '_' + param[2] + '_' + param[3] + '_' +param[4] + '_' + param[5] + '.RData'

        workflow.add_task_gridable('Rscript Load_Simulator.R' + ' --metadataType ' + param[6] + ' --nMicrobes '
        + param[0] + ' --nSamples ' + param[1] + ' --spikePerc ' + param[2] + ' --spikeStrength ' + param[3] + ' --nMetadata '
        + param[4] + ' --spikeMeta ' + param[5]  + ' --noZeroInflate' + ' --workingDirectory' + ' ' + str(working_directory)
        ,time=args.time,mem=args.memory,cores=args.threads,target=output)

    else:
        output = 'ZeroInflate' + '_' + param[6] + '_' + param[0] + '_' + param[1] + '_' + param[2] + '_' + param[3] + '_' +param[4] + '_' + param[5] + '.RData'

        workflow.add_task_gridable('Rscript Load_Simulator.R' + ' --metadataType ' + param[6] + ' --nMicrobes '
        + param[0] + ' --nSamples ' + param[1] + ' --spikePerc ' + param[2] + ' --spikeStrength ' + param[3] + ' --nMetadata '
        + param[4] + ' --spikeMeta ' + param[5] + ' --workingDirectory' + ' ' + str(working_directory)
        ,time=args.time,mem=args.memory,cores=args.threads,target=output)



workflow.go()
