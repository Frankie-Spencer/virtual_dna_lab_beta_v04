import os


def make_rnf_file(k_list, run_rnf_data, temper_data, run_time_data):

    data = []

    xml = '-xml {}\\{}.xml'.format(run_rnf_data['input_file_path_raw'], run_rnf_data['input_file_name'])
    utl = '-utl 1000'
    gdat = '-o {}.gdat'.format(run_rnf_data['input_file_name'])
    period = str(float(run_time_data[1]) / float(run_time_data[2])) if temper_data else run_time_data[1]
    dump = '-dump [{}:{}:{}]->{}'.format(run_time_data[0], period,
                                         run_time_data[3], run_rnf_data['dump_dir_path'])

    data.append('# written by VDNA ----- begin' + '\n')
    data.append('\n')

    data.append(xml + '\n')
    data.append(utl + '\n')
    data.append(gdat + '\n')
    data.append(dump + '\n')
    data.append('\n')

    data.append('begin\n')
    data.append('\n')

    for ii in k_list.items():
        data.append('set ' + ii[0] + ' ' + ii[1] + '\n')

    data.append('\n')
    data.append('update\n')

    if temper_data:
        ann_temp_data = [temper_data[0]]
        n_stages = int(temper_data[3])

        for stage_n in range(1, n_stages):
            if float(temper_data[0]) > float(temper_data[2]):
                new_temp = float(ann_temp_data[-1]) - float(temper_data[1])
                if new_temp >= float(temper_data[2]):
                    ann_temp_data.append(new_temp)
                else:
                    ann_temp_data.append(temper_data[2])

            elif float(temper_data[0]) < float(temper_data[2]):
                new_temp = float(ann_temp_data[-1]) + float(temper_data[1])
                if new_temp <= float(temper_data[2]):
                    ann_temp_data.append(new_temp)
                else:
                    ann_temp_data.append(temper_data[2])

        rts = []

        def save_data(e):
            rts.append('echo Temperature is Temp=' + str(e))
            rts.append('  set Temp ' + str(e))
            rts.append('  update')
            rts.append('  sim ' + str(run_time_data[1]) + ' ' + str(run_time_data[2]))
            rts.append('')

        for n_stage in ann_temp_data:
            save_data(n_stage)

        data.append('\n')

        for s in rts:
            new_line = str(s) + '\n'
            data.append(new_line)

        data.append('end\n')
        data.append('\n')
        data.append('# written by VDNA ----- end' + '\n')

    else:
        data.append('\n')
        data.append('sim ' + str(run_time_data[3]) + '\n')
        data.append('\n')
        data.append('end\n')
        data.append('\n')
        data.append('# written by VDNA ----- end' + '\n')

    f_name = 'current_simulation.rnf'
    with open('system_files/temp/' + f_name, 'w') as f:
        for item in data:
            f.write("%s" % item)
    f.close()

    rnf_file_path = os.path.abspath(__file__).rsplit('\\', 1)[0] + '\\temp\\' + f_name

    return rnf_file_path
