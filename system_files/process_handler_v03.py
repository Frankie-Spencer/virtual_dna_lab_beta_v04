import os
import shutil
import dump_to_species_converter_v25
import complexes_visualize_dev_v21_c_v05 as complex_visualize
from datetime import datetime

# destination directory folder names
save_species_file = 'species_files'
save_bngl = 'result_bngl'
save_visualize = 'visualize'


# prepare necessary prerequisite for main process
# this requires source and destination details as arguments
# this will create a destination folder named using input file name and current date-time
# recognise file format and convert if needed
# call complexes_visualize_dev_vXX by parsing the relevant arguments
def open_file(s, d, adv):

    def get_source(s):
        file_dir_name = s.rsplit('/', 1)
        file_name_full = file_dir_name[1]
        file_name = file_dir_name[1].rsplit('.', 1)[0]
        file_type = file_dir_name[1].rsplit('.', 1)[1]

        file = {'file_name_full': file_name_full,
                'file_name': file_name,
                'file_type': file_type
                }
        return file

    source = get_source(s)
    des_folder_name = datetime.now().strftime("%d-%m-%Y--%H.%M.%S")
    f_des_dir = source['file_name'] + '--' + '(' + des_folder_name + ')' + '/'
    des = d + '/' + f_des_dir if not d[-1] == '/' else d + f_des_dir

    try:
        os.makedirs(des, exist_ok=True)
        shutil.copytree('system_files/styles', des + 'styles')
    except:
        return ['Unable to process ✘, please try again!']

    if os.path.isdir(des):

        if s != '':

            try:
                if source['file_type'] == 'species':
                    source_path = shutil.copyfile(s, des + source['file_name'] + '_source_species' + '.species')
                    visualize = complex_visualize.complexes_vixualize(source_path, source['file_name'], des, adv)
                    if visualize:
                        complex_visualize.view(des + source['file_name'] + '_visualize.html')
                        return [True, source['file_name'], 'Species']

                elif source['file_type'] == '0':
                    source_path = des + source['file_name'] + '_source_species' + '.species'
                    dump_to_species_converter_v25.dump_to_species(s, source_path, source['file_name'], 'save_species')
                    visualize = complex_visualize.complexes_vixualize(source_path, source['file_name'], des, adv)
                    if visualize:
                        complex_visualize.view(des + source['file_name'] + '_visualize.html')
                        return [True, source['file_name'], 'Dump']

            except:
                shutil.rmtree(des)
                return ['Unable to process ✘, input file corrupted or not acceptable format!']

        else:
            return ['No input file given ✘']
    else:
        return ['Error occurred ✘, process aborted!']
