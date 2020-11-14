import anadama2
#Now create workflow
workflow = anadama2.Workflow()
workflow.add_argument("threads", desc="the number of threads for each task", default=4)
workflow.add_argument("time",desc="the time in minutes for each task",default=3*60)
workflow.add_argument("memory",desc="the memory in megabytes for each task",default=12*1024)
workflow.add_argument('parameters',desc="The filename of the parameters file",default='')
args = workflow.parse_args()



#import parameters from
paramscript = open(str(args.parameters),'r')
parameter_inputs = []
for line in paramscript:
    split_line = line.split()
    parameter_inputs.append(split_line)

method_list = parameter_inputs[0][1:]
metadata_type = parameter_inputs[1][1:]
n_microbes = parameter_inputs[2][1:]
n_samples = parameter_inputs[3][1:]
spike_perc = parameter_inputs[4][1:]
spike_strength = parameter_inputs[5][1:]
n_metadata = parameter_inputs[6][1:]
spike_meta = parameter_inputs[7][1:]
noZeroInflate = parameter_inputs[8][1:]
working_directory = parameter_inputs[9][1:]


#Create list of sets of parameters so we can check the output and resubmit failed jobs
param_list_evaluation = []
for microbe in n_microbes:
    for sample in n_samples:
        for perc in spike_perc:
            for size in spike_strength:
                for metadata in n_metadata:
                    for meta in spike_meta:
                        for metatype in metadata_type:
                            for inflate in noZeroInflate:
                                for method in method_list:
                                if metatype == 'UVA' or metatype == 'UVB' or metatype == 'UVC':

                                    param = (microbe,sample,perc,size,str(1),str(1),metatype,bool(inflate),method)
                                    param_list_generation.append(param)
                                else:
                                    param = (microbe,sample,perc,size,metadata,meta,metatype,bool(inflate),method)
                                    param_list_generation.append(param)

param_list_generation = []
for microbe in n_microbes:
    for sample in n_samples:
        for perc in spike_perc:
            for size in spike_strength:
                for metadata in n_metadata:
                    for meta in spike_meta:
                        for metatype in metadata_type:
                            for inflate in noZeroInflate:
                                if metatype == 'UVA' or metatype == 'UVB' or metatype == 'UVC':

                                    param = (microbe,sample,perc,size,str(1),str(1),metatype,bool(inflate))
                                    param_list.append(param)
                                else:
                                    param = (microbe,sample,perc,size,metadata,meta,metatype,bool(inflate))
                                    param_list_generation.append(param)


param_list_evaluation = list(set(param_list_evaluation))
param_list_generation = list(set(param_list_generation))
