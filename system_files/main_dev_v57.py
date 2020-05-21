from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QDialog
from sys_cache.cache import read_s_loc, write_s_loc, read_d_loc, write_d_loc, read_browser_loc
import os, time, threading, webbrowser, shutil
from dump_to_species_converter_v25 import dump_to_species
import process_handler_v03
from extract_ssdna_from_data_v02 import extract_ssdna
from make_rnf_v07 import make_rnf_file
from math import ceil
from multiprocessing.pool import ThreadPool


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        global write_edited_amount_g, write_edited_seq_g, advanced_criteria, advanced_criteria_text, ssdna_comp_list, \
               imp_seq_temp, type_dic_ssdna_comp, vis_adv_ranges_g, \
               run_kt_list_g, run_temper_data_g, run_time_data_g, run_annealing_g, run_kt_g, \
               write_kt_list_g, write_kt_status_g, \
               write_n_run_kt_validate, write_n_run_kt_validated_g

        run_kt_list_g, run_temper_data_g, run_time_data_g = {}, ['', '', ''], ['', '', '']
        run_kt_g, run_annealing_g = False, False

        write_kt_list_g, write_kt_status_g = {}, False

        write_n_run_kt_validate, write_n_run_kt_validated_g = 'wktv', ''

        type_dic_ssdna_comp = {1: 'ssDNA(s)    ', 2: 'complex(es) ', 'spaces': int}

        write_edited_amount_g, write_edited_seq_g = '', ''
        vis_adv_ranges_g = ''
        advanced_criteria, advanced_criteria_text = [], []
        imp_seq_temp = []
        ssdna_comp_list = []

        begin_species, end_species, ti = 'begin species', 'end species', [0]
        begin_parameters, end_parameters = 'begin parameters', 'end parameters'
        write_b_source = 'wbs'
        run_b_source, run_b_des = 'rbs', 'rbd'
        visual_b_source, visual_b_des = 'vbs', 'vbd'
        write_tab, run_tab, visual_tab = 'write', 'run', 'vis'
        write_the_seq, write_n_seq = 'wts', 'rns'
        run_start_time, run_n_dumps, run_sim_end, run_run = 'rst', 'rnd', 'rse', 'rr'
        import_seq = 'isf'
        species_type, dump_type, bngl_type = 'species', '0', 'bngl'
        write_edit_seq, write_edit_amount = 'res', 'wea'
        vis_adv_ranges = 'var'
        info_request = 'inr'

        def run_in_thread(fn):
            def run(*k, **kw):
                t = threading.Thread(target=fn, args=k, kwargs=kw)
                t.start()
                return t

            return run

        def get_user_consent(consent_for, s):
            consent_check = QMessageBox()

            con_dic = {'write_reset': 'Resetting will remove all created sequences. '
                                      'This action cannot be undone!\n' +
                                      '\n'
                                      'Click "Yes" to reset for a fresh start\n'
                                      'Click "No" to return back',

                       'delete_folder': 'Are you sure you want to permanently delete this folder?\n' +
                                        '\n' +
                                        s[0],

                       'delete_seq_comp': 'Are you sure you want to permanently delete this ' + s[1] + '?\n' +
                                          '\n' +
                                          s[0],

                       'simulate_line': 'Unable to process the given input file, because it contains the '
                                        'following syntax.\n' +
                                        '\n' +
                                        'line ' + s[1] + ':\n' +
                                        s[0] + '\n' +
                                        'Click "Yes" if you consent the program to comment out this line\n' +
                                        'Click "No" if you wish to do it manually'}

            if consent_for == 'write_reset' or consent_for == 'delete_folder' or 'delete_seq_comp':
                consent_check.setWindowTitle('Warning!')
                consent_check.setText(con_dic[consent_for])
                consent_check.setStandardButtons(QMessageBox.Yes | QMessageBox.No)  # | QtWidgets.QMessageBox.Cancel)
                consent = consent_check.exec_()

                return consent

            elif consent_for == 'simulate_line':
                self.label_messages.setText('Input file validation error ✘')

                consent_check.setWindowTitle('BNGL file validation error')
                consent_check.setText(con_dic[consent_for])
                consent_check.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                consent = consent_check.exec_()

                return consent

        def update_history():
            sele_folder_dir = self.lineEdit_visual_browse_des.text()

            if os.path.isdir(sele_folder_dir):
                self.listWidget_visual_history_list.clear()
                self.pushButton_visual_view_comp.setDisabled(True)
                self.pushButton_visual_view_bngl.setDisabled(True)
                self.pushButton_visual_view_source.setDisabled(True)
                self.pushButton_visual_delete.setDisabled(True)
                sele_folder = list(next(os.walk(sele_folder_dir))[:-1])

                c_time_dic_all = {}
                for dir in sele_folder[1]:
                    c_time_dic = {dir: os.stat(sele_folder[0] + '/' + dir).st_ctime}
                    c_time_dic_all.update(c_time_dic)

                sorted_all = dict(sorted(c_time_dic_all.items(), key=lambda kv: kv[1], reverse=True))

                for x in sorted_all:
                    if x != '':
                        self.listWidget_visual_history_list.insertItem(len(sorted_all), '✅  ' + str(x))
            else:
                self.listWidget_visual_history_list.clear()
                self.pushButton_visual_view_comp.setDisabled(True)
                self.pushButton_visual_view_bngl.setDisabled(True)
                self.pushButton_visual_view_source.setDisabled(True)
                self.pushButton_visual_delete.setDisabled(True)

        def get_set_file_name(a):
            source_path = self.lineEdit_write_browse_source.text()
            cus_file_name = self.lineEdit_custom_file_name.text()

            if os.path.isfile(source_path):
                r_source = source_path.rsplit('/', 1)
                source_name = r_source[1].rsplit('.', 1)[0]
                f_extension = r_source[1].rsplit('.', 1)[1]

                cus_checked = self.checkBox_custom_file_name.isChecked()

                default_file_name = source_name + '_VDNA' + '.' + f_extension
                c_file_name = cus_file_name + '.' + f_extension

                name_dic = {True: c_file_name,
                            False: default_file_name}

                rewrite = r_source[0] + '/' + name_dic[cus_checked]

                if a == info_request:
                    return [source_path, rewrite, name_dic[cus_checked]]
                else:
                    if not cus_checked:
                        self.lineEdit_custom_file_name.setText(default_file_name[:-5])
            else:
                self.lineEdit_custom_file_name.setText('')

        def input_file_validations(file_path, validate_call_from):
            global run_kt_list_g, write_kt_status_g, write_kt_list_g, write_kt_status_g

            source_file_dic = {write_b_source: self.lineEdit_write_browse_source.text(),
                               run_b_source: self.lineEdit_run_browse_source.text(),
                               import_seq: file_path,
                               run_run: file_path,
                               visual_b_source: self.lineEdit_visual_browse_source.text()}

            source_file_req = source_file_dic[validate_call_from]

            def validate_file(f_path, v_call_from):

                file_exists = os.path.isfile(f_path)

                if file_exists:
                    if v_call_from == write_b_source:
                        check_b_n_e = []
                        with open(f_path, 'r') as f:
                            a = f.readlines()
                            for i in a:
                                check_syn = i.rstrip('\n').rstrip()
                                if check_syn == begin_species:
                                    check_b_n_e.append(begin_species)

                                elif check_syn == end_species:
                                    check_b_n_e.append(end_species)

                        if begin_species in check_b_n_e and end_species in check_b_n_e:
                            return 'ready'
                        else:
                            none = check_b_n_e[0] if check_b_n_e != [] else 'none'
                            return none

                    elif v_call_from == run_b_source or v_call_from == run_run:
                        def re_write_input_file(c):

                            data = []

                            with open(c, 'r') as f:
                                a = f.readlines()
                                for i in a:
                                    if not i.startswith('simulate'):
                                        data.append(i)

                                    elif i.startswith('simulate'):
                                        new_line = '# ' + i
                                        data.append(new_line)
                            f.close()

                            with open(c, 'w') as f:
                                for item in data:
                                    f.write("%s" % item)
                            f.close()
                            return True

                        check_simulate = []
                        with open(f_path, 'r') as f:
                            a = f.readlines()
                            line_n = 0
                            for i in a:
                                line_n += 1
                                if i.startswith('simulate'):
                                    check_simulate.append(i)
                                    check_simulate.append(line_n)
                        f.close()

                        if check_simulate == []:
                            return True

                        elif check_simulate != []:
                            consent = get_user_consent('simulate_line', [check_simulate[0], str(check_simulate[1])])
                            if consent == 16384:
                                re_write_input_file(f_path)
                                return True
                            elif consent == 65536:
                                return False

                    elif v_call_from == import_seq or v_call_from == visual_b_source:
                        check_type = file_path.rsplit('/', 1)[1].rsplit('.', 1)[1]
                        imported_seqs = []
                        global run_status
                        run_status = None

                        @run_in_thread
                        def import_sep_thread(file_type):
                            global imp_seq_temp, run_status

                            def fr_dump():
                                dump_convert = dump_to_species(f_path, '', '', 'read_dump')
                                if dump_convert != False:
                                    get_ssdna = extract_ssdna('', 'read_species', dump_convert)
                                    return get_ssdna
                                else:
                                    return False

                            r_type_dic = {'0': fr_dump(),
                                          'species': extract_ssdna('', 'read_species', imported_seqs)}

                            import_run = list(r_type_dic[file_type])

                            if import_run != False:
                                imp_seq_temp = import_run
                                run_status = True
                            else:
                                run_status = False

                        with open(f_path, 'r') as f:
                            if check_type == species_type:
                                a = f.readlines()
                                for i in a:
                                    if i.startswith('N('):
                                        seq = str(i.rstrip('\n').rstrip()).rsplit('  ', 1)
                                        if seq[0].startswith('N(') and seq[1].isdigit():
                                            imported_seqs.append(seq[0] + '  ' + seq[1])
                                        else:
                                            return False

                                if v_call_from == import_seq:
                                    import_sep_thread('species')
                                    while run_status == None:
                                        time.sleep(0.1)

                                    if run_status == True:
                                        return True
                                    else:
                                        return False

                                elif v_call_from == visual_b_source:
                                    return [True, species_type]

                            elif check_type == dump_type:
                                if v_call_from == visual_b_source:
                                    return [True, dump_type]

                                elif v_call_from == import_seq:
                                    import_sep_thread('0')
                                    while run_status == None:
                                        time.sleep(0.1)

                                    if run_status == True:
                                        return True
                                    else:
                                        return False
                                else:
                                    return False

                else:
                    return 'file_not_exists'

            validate_check = validate_file(file_path, validate_call_from)

            if validate_call_from == write_b_source:

                if validate_check == 'ready':
                    write_s_loc(write_b_source, file_path)
                    get_set_file_name('')
                    write_kt_list_g = {}
                    write_kt_status_g = False
                    self.checkBox_write_advanced.setChecked(False)

                elif validate_check == end_species:
                    self.label_messages.setText(
                        'Source file error, "begin species" line not stated in source file ✘')

                elif validate_check == begin_species:
                    self.label_messages.setText('Source file error, "end species" line not stated in source file ✘')

                elif validate_check == 'none':
                    self.label_messages.setText('Source file error, '
                                                '"begin species" and "end species" '
                                                'lines not stated in source file ✘')

                elif not os.path.isfile(file_path) and os.path.isfile(source_file_req):
                    self.label_messages.setText('New source file not selected, returned to previous ✔')

                elif not os.path.isfile(file_path) and not os.path.isfile(source_file_req):
                    self.label_messages.setText('Source file not selected ✘')

            elif validate_call_from == run_b_source or validate_call_from == run_run:
                if validate_check == True:
                    write_s_loc(run_b_source, file_path)
                    self.label_messages.setText('BNGL file selected ✔')

                    if validate_call_from == run_b_source:
                        run_kt_list_g = {}
                        run_kt_status_g = False
                        self.checkBox_run_advanced.setChecked(False)

                elif validate_check == False:
                    self.lineEdit_run_browse_source.setText('')
                    self.label_messages.setText('BNGL file not selected ✘')

                elif not os.path.isfile(file_path) and os.path.isfile(source_file_req):
                    if validate_file(source_file_req, validate_call_from):
                        write_s_loc(run_b_source, source_file_req)
                        self.label_messages.setText('New BNGL file not selected, returned to previous ✔')

                    else:
                        self.lineEdit_run_browse_source.setText('')
                        self.label_messages.setText('Please edit the BNGL file and try again!')

                elif not os.path.isfile(file_path) and not os.path.isfile(source_file_req):
                    self.lineEdit_run_browse_source.setText('')
                    self.label_messages.setText('Source file not selected ✘')

            elif validate_call_from == import_seq:
                if validate_check == True:
                    write_s_loc(import_seq, file_path)

                else:
                    self.label_messages.setText('Species/Dump file error, sequence(s) import unsuccessful ✘')

            elif validate_call_from == visual_b_source:
                dic_type = {'species': 'Species', '0': 'Dump'}
                if validate_check[0] == True:
                    write_s_loc(visual_b_source, file_path)
                    self.label_messages.setText(
                        'Selected source file recognized as "' + dic_type[validate_check[1]] + '" file ✔')

                elif not os.path.isfile(file_path) and os.path.isfile(source_file_req):
                    self.label_messages.setText('New source file not selected, returned to previous ✔')

                elif not os.path.isfile(file_path) and not os.path.isfile(source_file_req):
                    self.label_messages.setText('Source file not selected ✘, please try again!')

            return validate_check

        def validate_and_update_buttons(tab):

            if tab == write_tab:
                created_list = self.listWidget_write_list_created
                list_all = [created_list.item(index).text() for index in range(created_list.count())]
                cus_file_name = self.lineEdit_custom_file_name.text()
                cus_checked = len(cus_file_name) > 0 if self.checkBox_custom_file_name.isChecked() else True
                input_validations = input_file_validations(self.lineEdit_write_browse_source.text(),
                                              write_b_source) == 'ready'

                if len(list_all) > 0:
                    self.pushButton_write_reset_all.setEnabled(True)
                    if input_validations and cus_checked:
                        self.pushButton_write_submit_with.setEnabled(True)
                        self.pushButton_write_submit_without.setEnabled(True)
                    else:
                        self.pushButton_write_submit_with.setDisabled(True)
                        self.pushButton_write_submit_without.setDisabled(True)
                else:
                    self.pushButton_write_reset_all.setDisabled(True)
                    self.pushButton_write_submit_with.setDisabled(True)
                    self.pushButton_write_submit_without.setDisabled(True)

                if created_list.currentRow() == -1:
                    self.pushButton_write_delete.setDisabled(True)
                    self.pushButton_write_edit.setDisabled(True)
                else:
                    self.pushButton_write_delete.setEnabled(True)
                    self.pushButton_write_edit.setEnabled(True)

                if len(self.lineEdit_write_sequence.text()) >= 2 and len(
                        self.lineEdit_write_amount.text()) >= 1:
                    self.pushButton_write_add_sequence.setEnabled(True)
                else:
                    self.pushButton_write_add_sequence.setDisabled(True)

                if input_validations:
                    self.pushButton_write_advanced.setEnabled(True)
                    self.checkBox_write_advanced.setDisabled(False)
                else:
                    self.pushButton_write_advanced.setDisabled(True)
                    self.checkBox_write_advanced.setDisabled(True)

            elif tab == run_tab:
                start_time = self.lineEdit_run_start_time.text()
                end_time = self.lineEdit_run_sim_end.text()
                n_dumps = self.lineEdit_run_n_dumps.text()
                source_file_path = self.lineEdit_run_browse_source.text()
                des_dir_path = self.lineEdit_run_browse_des.text()

                st_time = str(start_time)[-1] if len(str(start_time)) > 0 else None
                ed_time = str(end_time)[-1] if len(str(end_time)) > 0 else None

                if os.path.isfile(source_file_path):
                    if source_file_path.rsplit('/', 1)[1].rsplit('.', 1)[1] == bngl_type:
                        self.pushButton_run_advanced.setEnabled(True)
                        self.checkBox_run_advanced.setDisabled(False)
                        if os.path.isdir(des_dir_path):
                            if start_time != '' and end_time != '' and n_dumps != '':
                                if not st_time.endswith('.') and not ed_time.endswith('.'):
                                    if float(start_time) != float(end_time):
                                        if float(start_time) < float(end_time):
                                            if float(end_time) > 0:
                                                if int(n_dumps) > 0:
                                                    self.pushButton_run_run.setEnabled(True)
                                                    self.label_messages.setText('Program ready to run!')
                                                    return True

                                                else:
                                                    self.pushButton_run_run.setDisabled(True)
                                                    self.label_messages.setText('"# dump files" '
                                                                                'must be a positive value ✘')

                                            else:
                                                self.pushButton_run_run.setDisabled(True)
                                                self.label_messages.setText('"End time" must be a positive value ✘')

                                        else:
                                            self.pushButton_run_run.setDisabled(True)
                                            self.label_messages.setText('"Start time" must be less than "End time" ✘')

                                    else:
                                        self.pushButton_run_run.setDisabled(True)
                                        self.label_messages.setText('"Start time" and "End time" cannot be same ✘')

                                else:
                                    self.pushButton_run_run.setDisabled(True)
                                    self.label_messages.setText('Error found on some fields, please check ✘')

                            else:
                                self.pushButton_run_advanced.setEnabled(True)
                                self.pushButton_run_run.setDisabled(True)
                                self.label_messages.setText('Any of the below fields cannot be empty ✘')

                        else:
                            self.pushButton_run_advanced.setEnabled(True)
                            self.pushButton_run_run.setDisabled(True)
                            self.label_messages.setText('Save simulation outputs directory not selected ✘')

                    else:
                        self.checkBox_run_advanced.setDisabled(True)
                        self.pushButton_run_advanced.setDisabled(True)
                        self.pushButton_run_run.setDisabled(True)
                        self.label_messages.setText('Input file error ✘, only "bngl" files accepted!')

                else:
                    self.checkBox_run_advanced.setDisabled(True)
                    self.pushButton_run_run.setDisabled(True)
                    self.pushButton_run_advanced.setDisabled(True)
                    self.label_messages.setText('Input BNGL file is not selected ✘')

            elif tab == visual_tab:
                created_list = self.listWidget_visual_history_list
                source_file_path = self.lineEdit_visual_browse_source.text()
                des_dir_path = self.lineEdit_visual_browse_des.text()

                accepted_types = [dump_type, species_type]

                if created_list.currentRow() == -1:
                    self.pushButton_visual_view_comp.setDisabled(True)
                    self.pushButton_visual_view_bngl.setDisabled(True)
                    self.pushButton_visual_delete.setDisabled(True)
                else:
                    if os.path.isdir(des_dir_path):
                        history_list = self.listWidget_visual_history_list
                        list_all = [history_list.item(index).text() for index in range(history_list.count())]
                        selection = list_all[history_list.currentRow()][3:]
                        selection_link = self.lineEdit_visual_browse_des.text() + '/' + selection

                        b_status = [False, False, False]
                        for file in os.listdir(selection_link):

                            if file.endswith('.html'):
                                b_status[0] = True

                            if file.endswith('.species') and '_result_species' in file:
                                b_status[1] = True

                            if file.endswith('.species') and '_source_species' in file:
                                b_status[2] = True

                        self.pushButton_visual_view_comp.setEnabled(b_status[0])
                        self.pushButton_visual_view_bngl.setEnabled(b_status[1])
                        self.pushButton_visual_view_source.setEnabled(b_status[2])

                        self.pushButton_visual_delete.setEnabled(True)

                if not os.path.isfile(source_file_path):
                    self.pushButton_visual_visualize.setDisabled(True)
                    self.label_messages.setText('Input dump/species file is not selected ✘')

                elif os.path.isfile(source_file_path) and \
                        source_file_path.rsplit('/', 1)[1].rsplit('.', 1)[1] not in accepted_types:
                    self.pushButton_visual_visualize.setDisabled(True)
                    self.label_messages.setText('Input file error ✘, only ".speceis" or ".0" files accepted!')

                    '''
                    elif input_file_type in accepted_types:
                        self.pushButton_visual_visualize.setEnabled(True)
                        self.label_messages.setText(
                            'Selected source file recognized as "' + dic_type[input_file_type] + '" file ✔')
                    '''
                elif not os.path.isdir(des_dir_path):
                    self.pushButton_visual_visualize.setDisabled(True)
                    self.label_messages.setText('Save simulation outputs directory not selected ✘')

                else:
                    self.pushButton_visual_visualize.setEnabled(True)
                    self.label_messages.setText('Program ready to run!')
                    return True

        def on_edit_source_path(call_from, tab):
            source_file_dic = {write_b_source: self.lineEdit_write_browse_source.text(),
                               run_b_source: self.lineEdit_run_browse_source.text(),
                               visual_b_source: self.lineEdit_visual_browse_source.text()}

            file_path = source_file_dic[call_from]

            if os.path.isfile(file_path):
                input_file_validations(file_path, call_from)

            validate_and_update_buttons(tab)

        def browse_path(browse_call_from, tab):
            source_file_dic = {write_b_source: self.lineEdit_write_browse_source,
                               run_b_source: self.lineEdit_run_browse_source,
                               visual_b_source: self.lineEdit_visual_browse_source}

            browse_sources = [write_b_source, run_b_source, visual_b_source, import_seq]

            browse_des = {run_b_des: [self.lineEdit_run_browse_des, 'Save simulation outputs selected ✔'],
                          visual_b_des: [self.lineEdit_visual_browse_des, 'Save analysis outputs selected ✔']}

            win_heading_and_file_type = {write_b_source: ['Browse source BNGL file', 'BNGL(*.bngl)'],
                                         run_b_source: ['Browse source BNGL file', 'BNGL(*.bngl)'],
                                         visual_b_source: ['Browse Dump/Species file', 'Dump/Species(*.0 *.species)'],
                                         visual_b_des: 'Save analysis outputs',
                                         run_b_des: 'Save simulation outputs',
                                         import_seq: ['Import sequences from Dump/Species file',
                                                      'Dump/Species(*.0 *.species)']}

            def validate_path(link):
                chars = ['#']
                chars_in_link = []

                for c in link:
                    if c in chars:
                        chars_in_link.append(c)
                if chars_in_link == []:
                    return [True, '']

                elif chars_in_link != []:
                    chars_list = '   ' + '  '.join(chars_in_link) + '   '
                    return [False, chars_list]

            if browse_call_from in browse_sources:

                path = str(os.path.dirname(str(read_s_loc(browse_call_from))))
                file_path, _ = QFileDialog.getOpenFileName(None, win_heading_and_file_type[browse_call_from][0],
                                                           path, win_heading_and_file_type[browse_call_from][1])

                if len(file_path) > 0:
                    validate_path_status = validate_path(file_path)

                    if validate_path_status[0]:
                        if not browse_call_from == import_seq:
                            source_file_dic[browse_call_from].setText(file_path)

                        elif browse_call_from == import_seq:
                            validate_status = input_file_validations(file_path, browse_call_from)
                            return validate_status

                        validate_and_update_buttons(tab)

                    else:
                        self.label_messages.setText('Path cannot include special character(s) ' + validate_path_status[1])

            elif browse_call_from in browse_des:
                path = str(read_d_loc(browse_call_from))
                dir_path = QFileDialog.getExistingDirectory(None, win_heading_and_file_type[browse_call_from], path)

                cur_des_path = self.lineEdit_visual_browse_des.text()

                validate_path_status = validate_path(dir_path)

                if validate_path_status[0]:
                    if os.path.isdir(dir_path):
                        browse_des[browse_call_from][0].setText(dir_path)
                        write_d_loc(browse_call_from, dir_path)

                    elif not os.path.isdir(dir_path) and os.path.isdir(cur_des_path):
                        self.label_messages.setText('New destination folder not selected, returned to previous ✔')

                    elif not os.path.isdir(dir_path) and not os.path.isdir(cur_des_path):
                        self.label_messages.setText('Destination folder not selected ✘, please try again!')

                    update_history() if browse_call_from == visual_b_des else None
                    self.label_messages.setText(browse_des[browse_call_from][1])

                    validate_and_update_buttons(tab)

                else:
                    self.label_messages.setText('Path cannot include special character(s) ' + validate_path_status[1])

        @run_in_thread
        def run_visualisation(a):
            s_file_path = self.lineEdit_visual_browse_source.text()
            d_dir_path = self.lineEdit_visual_browse_des.text()
            get_adv = advanced_criteria if len(
                advanced_criteria) != 0 and self.checkBox_visual_advanced.isChecked() else None
            run = browse_and_parse_v03.open_file(s_file_path, d_dir_path, get_adv)

            if isinstance(run[0], bool):
                self.label_messages.setText(run[2] + ' file --- ' + run[1] + ' --- processed successfully ✔')
                update_history()
            elif isinstance(run[0], str):
                self.label_messages.setText(run[0])
            else:
                return None

        def pop_view(a, b):
            history_list = self.listWidget_visual_history_list
            list_all = [history_list.item(index).text() for index in range(history_list.count())]
            selection = list_all[history_list.currentRow()][3:]
            selection_link = self.lineEdit_visual_browse_des.text() + '/' + selection

            try:
                for file in os.listdir(selection_link):
                    if file.endswith(a):
                        if b == 'html':
                            f_path = selection_link + '/' + file
                            chrome_path = read_browser_loc() + ' %s'
                            x = lambda: webbrowser.get(chrome_path).open_new(f_path)
                            t = threading.Thread(target=x)
                            t.start()

                        elif b == 's_sp':
                            if '_source_species' in file:
                                f_path = selection_link + '/' + file
                                os.startfile(f_path)

                        elif b == 'r_sp':
                            if '_result_species' in file:
                                f_path = selection_link + '/' + file
                                os.startfile(f_path)
            except:
                self.label_messages.setText('Some error occurred, please try again!')

        def delete_folder():
            history_list = self.listWidget_visual_history_list
            list_all = [history_list.item(index).text() for index in range(history_list.count())]
            selection = list_all[history_list.currentRow()][3:]
            selection_link = self.lineEdit_visual_browse_des.text() + '/' + selection

            folder_name = selection_link.rsplit('/', 1)[1]

            if os.path.isdir(selection_link):
                if os.path.isdir(selection_link):
                    warning = get_user_consent('delete_folder', [folder_name, ''])

                    if warning == 16384:
                        shutil.rmtree(selection_link)
                        self.label_messages.setText(folder_name + ' --- deleted ✔')
                        update_history()

                    else:
                        self.label_messages.setText('No changes made!')
            else:
                self.label_messages.setText('This folder does not exists anymore!')

        def create_bngl_sequence(ent):
            dna_seq = ['A', 'T', 'G', 'C']

            seq = ''.join(c.upper() for c in ent if c.upper() in dna_seq) if len(ent) != 0 else None

            if seq:
                bngl_syn = 'N(b~' + seq[0] + ',5,3!1,W)'

                if len(seq) > 2:
                    for i in range(1, len(seq) - 1):
                        bngl_syn += '.N(b~{},5!{},3!{},W)'.format(seq[i], str(i), str(i + 1))

                    bngl_syn += '.N(b~{},5!{},3,W)'.format(seq[-1], str(len(seq) - 1))

                elif len(seq) == 2:
                    bngl_syn += '.N(b~{},5!{},3,W)'.format(seq[-1], str(1))

                else:
                    self.label_messages.setText('Minimum two nucleotides required!')
                    return None

                return bngl_syn

            else:
                self.label_messages.setText('Minimum two nucleotides required!')

            validate_and_update_buttons(write_tab)

        def update_created_list():
            global ssdna_comp_list, type_dic_ssdna_comp

            if ssdna_comp_list != []:
                mx_size = len(str(max([k[1] for k in ssdna_comp_list])))
                type_dic_ssdna_comp['spaces'] = mx_size
                self.listWidget_write_list_created.clear()

                for seqs in reversed(ssdna_comp_list):
                    s_text = '{} || {}{} ||  {}'.format(type_dic_ssdna_comp[seqs[0]], seqs[1],
                                                        ' ' * int(mx_size - len(str(seqs[1]))),
                                                        seqs[2])
                    self.listWidget_write_list_created.insertItem(0, s_text)
            else:
                type_dic_ssdna_comp['spaces'] = int
                self.listWidget_write_list_created.clear()

            validate_and_update_buttons(write_tab)

        def add_sequence():
            global ssdna_comp_list

            seq = self.lineEdit_write_sequence.text()
            am = self.lineEdit_write_amount.text() or '1'

            exist_items = [e[2] for e in ssdna_comp_list]

            if seq not in exist_items:
                seq_syn = create_bngl_sequence(seq)
                ll = [1, int(am), seq, seq_syn]
                ssdna_comp_list.insert(0, ll)

            elif seq in exist_items:
                item_index = exist_items.index(seq)
                ssdna_comp_list[item_index][1] += int(am)

            update_created_list()

        def import_sequences(call_from, tab):
            global imp_seq_temp, ssdna_comp_list

            browse_get_status = browse_path(call_from, tab)

            if browse_get_status == True:
                exist_items = ['|' + e[3] + '|' for e in ssdna_comp_list]

                for ll in reversed(imp_seq_temp):
                    u_item = '|' + ll[3] + '|'
                    if u_item not in exist_items:
                        ssdna_comp_list.insert(0, ll)
                        exist_items.insert(0, u_item)

                    elif u_item in exist_items:
                        item_index = exist_items.index(u_item)
                        ssdna_comp_list[item_index][1] += ll[1]

                self.label_messages.setText('Same ssDNA(s) or Complex(es) already exists, their amounts aggregated ✔')
                imp_seq_temp = []

                update_created_list()

        def validate_custom_name():
            custom_name = self.lineEdit_custom_file_name.text()
            not_allowed_chars = ['\\', '/', ':', '*', '?', '"', '<', '>', '|']
            n_a_char_msg = '\ / : * ? " < > |'
            new_name = ''

            for c in custom_name:
                if c not in not_allowed_chars:
                    new_name += c

            if custom_name != new_name:
                self.lineEdit_custom_file_name.setText(new_name)
                self.label_messages.setText('Custom name error,   ' + n_a_char_msg + '   characters not allowed ✘')
            else:
                self.label_messages.setText('')

            validate_and_update_buttons(write_tab)

        def validate_char(call_from):
            seq_in_var = self.lineEdit_write_sequence
            seq_in = {write_the_seq: seq_in_var.text(), write_edit_seq: write_edited_seq_g}
            cur_seq = seq_in[call_from]
            allowed = ('A', 'T', 'C', 'G')

            def set_seq_values(v):
                global write_edited_seq_g

                if cur_seq != v:
                    if call_from == write_the_seq:
                        seq_in_var.setText(v)
                    elif call_from == write_edit_seq:
                        write_edited_seq_g = v

            new_seq = ''
            for s in cur_seq:
                if len(cur_seq) != 0:
                    if s.upper() in allowed:
                        new_seq += s.upper()

            set_seq_values(new_seq)

            if not call_from == write_edit_seq:
                validate_and_update_buttons(write_tab)

        def validate_num(call_from, tab):

            place_dic = {write_n_seq: self.lineEdit_write_amount,
                         run_start_time: self.lineEdit_run_start_time,
                         run_n_dumps: self.lineEdit_run_n_dumps,
                         run_sim_end: self.lineEdit_run_sim_end,
                         write_edit_amount: write_edited_amount_g,
                         vis_adv_ranges: vis_adv_ranges_g,
                         write_n_run_kt_validate: write_n_run_kt_validated_g}

            child_window = [write_edit_amount, vis_adv_ranges, write_n_run_kt_validate]

            amount_in = str(place_dic[call_from].text()) if call_from not in child_window \
                else str(place_dic[call_from])

            only_int_type = [write_n_seq, run_n_dumps, write_edit_amount, vis_adv_ranges]
            only_float_type = [run_start_time, run_sim_end, write_n_run_kt_validate]

            if call_from in only_int_type:

                def set_amount_values(v):
                    global write_edited_amount_g, vis_adv_ranges_g

                    if amount_in != v:
                        if call_from not in child_window:
                            place_dic[call_from].setText(v)
                        elif call_from in child_window:
                            if call_from == write_edit_amount:
                                write_edited_amount_g = v
                            elif call_from == vis_adv_ranges:
                                vis_adv_ranges_g = v

                new_seq = ''
                if len(amount_in) >= 1:
                    for s in amount_in:
                        if s.isdigit():
                            new_seq += s

                    if len(new_seq) >= 1:
                        if int(new_seq) >= 1:
                            set_amount_values(str(int(new_seq)))
                        elif int(new_seq) == 0:
                            set_amount_values('1')
                    else:
                        set_amount_values('')

            if call_from in only_float_type:

                def set_amount_values(v):
                    global write_n_run_kt_validated_g

                    if amount_in != v:
                        if call_from not in child_window:
                            place_dic[call_from].setText(v)
                        elif call_from in child_window:
                            write_n_run_kt_validated_g = v

                new_seq_num = ''
                if len(amount_in) >= 1:
                    for s in amount_in:
                        if s.isdigit() or s == '.':
                            new_seq_num += s

                    coun_dots = new_seq_num.count('.')
                    new_seq_clean = ''
                    for ss in new_seq_num:
                        if ss != '.':
                            new_seq_clean += ss
                        elif ss == '.':
                            if '.' not in new_seq_clean:
                                if not coun_dots > 1:
                                    new_seq_clean += ss
                                elif coun_dots > 1:
                                    coun_dots -= 1

                    if len(new_seq_clean) == 1 and new_seq_clean == '.':
                        set_amount_values('')
                    else:
                        try:
                            new_seq_clean = int(new_seq_clean)
                        except:
                            pass

                        if not type(new_seq_clean) == int:
                            try:
                                if not new_seq_clean.startswith('.') and not new_seq_clean.endswith('.'):
                                    new_seq_clean_new = new_seq_clean.split('.')
                                    new_seq_clean = '{}.{}'.format(int(new_seq_clean_new[0]), new_seq_clean_new[1])
                                elif new_seq_clean.endswith('.'):
                                    new_seq_clean = str(int(new_seq_clean[:-1])) + '.'

                                elif new_seq_clean.startswith('.'):
                                    new_seq_clean = str(int(new_seq_clean[1:]))
                            except:
                                new_seq_clean = ''

                        set_amount_values(str(new_seq_clean))

            if call_from not in child_window:
                validate_and_update_buttons(tab)

        def write_bngl(a):
            source_info = get_set_file_name(info_request)
            source_path = source_info[0]
            ori_dic = {1: ", 5' - 3' orientation", 2: ''}

            if os.path.isfile(source_path):

                list_all = ssdna_comp_list

                if list_all != []:
                    rewrite_name = source_info[2]
                    rewrite = source_info[1]

                    data, q = [], []
                    type = {'with': ' with existing ssDNA sequence(s)', 'without': ' without existing sequence(s)'}
                    with open(source_path, 'r') as f:
                        l = f.readlines()
                        for i in l:
                            check_syn = i.rstrip('\n').rstrip()

                            if self.checkBox_write_advanced.isChecked() and end_parameters not in q:
                                if begin_parameters not in q:
                                    if check_syn != begin_parameters:
                                        data.append(i)

                                    elif check_syn == begin_parameters:
                                        data.append(check_syn + '\n')
                                        q.append(check_syn)
                                        data.append('\n')
                                        data.append('# written by VDNA ----- begin parameters' + '\n')
                                        data.append('\n')

                                        for s in write_kt_list_g.items():
                                            new_p_line = '  {} {}\n'.format(s[0], s[1])
                                            data.append(new_p_line)

                                        data.append('\n')
                                        data.append('# written by VDNA ----- end parameters' + '\n')
                                        data.append('\n')
                                        q.append('done_parameters')

                                elif begin_parameters in q and 'done_parameters' in q and end_parameters not in q:
                                    if check_syn != end_parameters:
                                        continue

                                    elif check_syn == end_parameters:
                                        data.append(check_syn + '\n')
                                        q.append(check_syn)

                            elif not self.checkBox_write_advanced.isChecked() or end_parameters in q:
                                if begin_species not in q:
                                    if check_syn != begin_species:
                                        data.append(i)

                                    elif check_syn == begin_species:
                                        data.append(check_syn + '\n')
                                        q.append(check_syn)
                                        data.append('\n')
                                        data.append('# written by VDNA ----- begin species' + '\n')
                                        data.append('\n')

                                        for s in list_all:
                                            spaces_amount = ' ' * int(type_dic_ssdna_comp['spaces'] - len(str(s[1])))
                                            '# {} || {}{} ||  {}{}\n'
                                            comment_comp_info = '# {} || {}{} ||  {}{}\n'.format(type_dic_ssdna_comp[s[0]],
                                                                                                 s[1], spaces_amount, s[2],
                                                                                                 ori_dic[s[0]])
                                            new_line = '{}  {}\n'.format(s[3], s[1])
                                            data.append(comment_comp_info)
                                            data.append(new_line)
                                            data.append('\n')
                                        data.append('# written by VDNA ----- end species' + '\n')
                                        data.append('\n')
                                        q.append('done_species')

                                elif begin_species in q and 'done_species' in q and end_species not in q:
                                    if check_syn != end_species and a == 'without':
                                        continue

                                    elif check_syn != end_species and a == 'with':
                                        data.append(i)

                                    elif check_syn == end_species:
                                        data.append(check_syn + '\n')
                                        q.append(check_syn)

                                elif begin_species in q and 'done_species' in q and end_species in q:
                                    data.append(i)
                    f.close()

                    with open(rewrite, 'w') as f:
                        for item in data:
                            f.write("%s" % item)
                    f.close()

                    self.label_messages.setText('BNGL file successfully saved as ' +
                                                '"' + rewrite_name + '"' +
                                                str(type[a]) + ' ✔')
                    self.lineEdit_run_browse_source.textChanged.disconnect()
                    self.lineEdit_run_browse_source.setText(rewrite)
                    write_s_loc(run_b_source, rewrite)
                    self.lineEdit_run_browse_source.textChanged.connect(lambda: validate_and_update_buttons(run_tab))
                else:
                    self.label_messages.setText('Please create ssDNA sequence(s) to write!')
            else:
                self.label_messages.setText('Please select the source file!')

        def delete_created_sequences(c_row):
            type_dic = {1: 'sequence', 2: 'complex'}
            seq = ssdna_comp_list[c_row]

            consent = get_user_consent('delete_seq_comp', [str(seq[2][:40] + ' ...'), type_dic[seq[0]]])

            if consent == 16384:
                del ssdna_comp_list[c_row]
                update_created_list()

            validate_and_update_buttons(write_tab)

        def reset_created_sequences_list():
            global ssdna_comp_list

            if get_user_consent('write_reset', ['', '']) == 16384:
                self.listWidget_write_list_created.clear()
                ssdna_comp_list = []
                self.label_messages.setText('List reset done ✔')

            validate_and_update_buttons(write_tab)

        def custom_name_set():
            checkbox_cfm_var = self.checkBox_custom_file_name
            if checkbox_cfm_var.isChecked():
                self.lineEdit_custom_file_name.setEnabled(True)
            else:
                self.lineEdit_custom_file_name.setDisabled(True)

            get_set_file_name('')

        def run_bngl():
            if input_file_validations(self.lineEdit_run_browse_source.text(), run_run) == True:
                start_time = self.lineEdit_run_start_time.text()
                end_time = self.lineEdit_run_sim_end.text()
                n_dumps = self.lineEdit_run_n_dumps.text()

                dump = {'start': start_time,
                        'period': str(round((float(end_time) / float(n_dumps)), 4)),
                        'end': end_time,
                        'steps': n_dumps}

                def get_source(s):
                    file_dir_name = s.rsplit('/', 1)
                    file_name_full = file_dir_name[1]
                    file_name = file_dir_name[1].rsplit('.', 1)[0]
                    file_type = file_dir_name[1].rsplit('.', 1)[1]
                    file_name_xml = file_name + '.xml'

                    file = {'file_dir': file_dir_name[0],
                            'file_name_full': file_name_full,
                            'file_name': file_name,
                            'file_type': file_type,
                            'file_name_xml': file_name_xml
                            }
                    return file

                input_info = get_source(self.lineEdit_run_browse_source.text())

                def creat_dump_folder():
                    dump_dir = self.lineEdit_run_browse_des.text()

                    def creat_folder(name):
                        nums = []
                        try:
                            existing = list(next(os.walk(dump_dir)))[1]
                            for i in existing:
                                n_num = i.rsplit('--', 1)
                                if name == n_num[0]:
                                    nums.append(int(n_num[1]))
                        except:
                            pass

                        new_name = name + '--' + (str(max(nums) + 1) if nums != [] else '1')
                        return new_name

                    dump_folder_name = creat_folder(input_info['file_name'])
                    dump_dir = dump_dir.replace('/', '\\') + '\\' + dump_folder_name

                    return {'dump_dir': dump_dir, 'dump_folder_name': dump_folder_name}

                dump_info = creat_dump_folder()

                def convert_link_address(d):
                    chars = ['\\', '/', '_', ':', '-', '.']
                    dir = ''

                    for c in d:
                        if c.isalpha() or c.isdigit() or c in chars:
                            dir += c
                        else:
                            dir += '^' + c
                    return dir

                local_path = os.path.abspath(__file__).rsplit('\\', 2)[0]
                input_file_path_raw = input_info['file_dir'].replace('/', '\\')
                dump_dir_path_converted = convert_link_address(dump_info['dump_dir'])
                dump_dir_path_for_rnf = dump_info['dump_dir'] + '/'

                nfsim_pl_path_converted = convert_link_address(local_path + '\\NFsim_v1.11\\bng2.pl')
                nfsim_exe_path_converted = convert_link_address(local_path + '\\NFsim_v1.11\\bin\\NFsim_MSWin32.exe')
                perl_path_converted = convert_link_address(local_path + '\\Perl64\\bin\\perl.exe')

                input_file_bngl = '{}\\{}'.format(input_file_path_raw, input_info['file_name_full'])
                input_file_xml = '{}\\{}'.format(input_file_path_raw, input_info['file_name_xml'])
                run_period = '[{}:{}:{}]'.format(dump['start'], dump['period'], dump['end'])

                com_change_dir = '"CD "{}\\"'.format(convert_link_address(input_file_path_raw))
                com_make_xml = '"{}" "{}" -xml "{}"'.format(perl_path_converted,
                                                            nfsim_pl_path_converted,
                                                            convert_link_address(input_file_bngl))
                com_make_dump_dir = 'md "{}\\"'.format(dump_dir_path_converted)
                com_run_nfsim_basic = '"{}" -utl 1000 -xml "{}" -dump "{}-^>{}/" -oSteps {} -sim {}"'.format(
                                nfsim_exe_path_converted,
                                convert_link_address(input_file_xml),
                                run_period,
                                dump_dir_path_converted,
                                dump['steps'],
                                end_time)

                main_run_command = 'START CMD /K {} && {} && {}'.format(com_change_dir,
                                                                        com_make_xml,
                                                                        com_make_dump_dir)

                run_rnf_data_dic = {'input_file_name': input_info['file_name'],
                                    'input_file_path_raw': input_file_path_raw,
                                    'dump_dir_path': dump_dir_path_for_rnf}

                def run_sim_on_thread(cmd_com):
                    x = lambda: os.system(cmd_com)
                    t = threading.Thread(target=x)
                    t.start()
                    self.label_messages.setText('Processing data on CMD!')

                if self.checkBox_run_advanced.isChecked() and run_kt_g and run_annealing_g:
                    make_rnf = make_rnf_file(run_kt_list_g, run_rnf_data_dic, run_temper_data_g, run_time_data_g)
                    com_run_rnf = '"{}" -rnf "{}"'.format(nfsim_exe_path_converted, convert_link_address(make_rnf))
                    cmd_command = '{} && {}'.format(main_run_command,
                                                    com_run_rnf)

                    run_sim_on_thread(cmd_command)

                elif self.checkBox_run_advanced.isChecked() and run_kt_g and not run_annealing_g:
                    run_time_data_basic = [dump['start'], dump['period'], '', dump['end']]
                    make_rnf = make_rnf_file(run_kt_list_g, run_rnf_data_dic, False, run_time_data_basic)
                    com_run_rnf = '"{}" -rnf "{}"'.format(nfsim_exe_path_converted, convert_link_address(make_rnf))
                    cmd_command = '{} && {}'.format(main_run_command,
                                                    com_run_rnf)

                    run_sim_on_thread(cmd_command)

                else:
                    cmd_command = '{} && {}'.format(main_run_command,
                                                    com_run_nfsim_basic)

                    run_sim_on_thread(cmd_command)

        def advanced_options():
            global p_o_n

            p_o_n = 0

            def get_add_criteria():
                range_from, range_to, unbound_from, unbound_to, criteria_list = lineEdit_range_from, \
                                                                                lineEdit_range_to, \
                                                                                lineEdit_unbound_from, \
                                                                                lineEdit_unbound_to, \
                                                                                listWidget_list_criteria
                l_dummy = [None, None, None, None, None, False]
                l_display = ''
                p_o_n_dic = {0: ['%', '0%', '100%', '%'], 1: ['n', 'Min.', 'Max.', '']}
                amount_from_item = range_from.text()
                amount_to_item = range_to.text()
                unbound_from_item = unbound_from.text()
                unbound_to_item = unbound_to.text()
                # p_o_n = p_or_n.text()
                range_text = 'ssDNA(s) range, from {} to {}'.format(
                    amount_from_item if amount_from_item else 'Min.',
                    amount_to_item if amount_to_item else 'Max.')

                if amount_from_item != '' or amount_to_item != '':
                    l_dummy[0] = amount_from_item
                    l_dummy[1] = amount_to_item

                if unbound_from_item != '' or unbound_to_item != '':
                    pon = p_o_n_dic[p_o_n]
                    l_dummy[2] = unbound_from_item
                    l_dummy[3] = unbound_to_item
                    l_dummy[4] = pon[0]
                    l_d_p = 'Unbound {}, from {} to {}'.format(pon[0], unbound_from_item + pon[3]
                    if unbound_from_item != '' else pon[1],
                                                               unbound_to_item + pon[3]
                                                               if unbound_to_item != '' else pon[2])
                    l_display = l_d_p

                new_item = '✅  {}{}'.format(range_text, '  |  ' + l_display if l_display != '' else '')

                return [l_dummy, new_item]

            def update_add_button():

                item = get_add_criteria()
                list_var = listWidget_list_criteria
                existing = [list_var.item(index).text() for index in range(list_var.count())]

                if item[0][0] or item[0][1]:
                    current_range_from = int(item[0][0]) if item[0][0] else 1
                    current_range_to = int(item[0][1]) if item[0][1] else float('inf')

                    if current_range_from <= current_range_to:
                        check_if_exists = item[1] in [i[:-14] if i.endswith('Highlighted') else i for i in existing]
                        if not item[0][2] and not item[0][3]:
                            if not check_if_exists:
                                pushButton_add_parameter.setEnabled(True)
                            else:
                                pushButton_add_parameter.setDisabled(True)

                        elif item[0][2] or item[0][3]:
                            unb_from_n = int(item[0][2]) if item[0][2] else 0
                            unb_to_n = int(item[0][3]) if item[0][3] else 100 if item[0][4] == '%' else float('inf')
                            if unb_from_n <= unb_to_n:
                                if not check_if_exists:
                                    pushButton_add_parameter.setEnabled(True)
                                else:
                                    pushButton_add_parameter.setDisabled(True)
                            else:
                                pushButton_add_parameter.setDisabled(True)
                        else:
                            pushButton_add_parameter.setDisabled(True)
                    else:
                        pushButton_add_parameter.setDisabled(True)
                else:
                    pushButton_add_parameter.setDisabled(True)

            def validate_filter_input(place_id):

                place_dic = {'range_from': lineEdit_range_from,
                             'range_to': lineEdit_range_to,
                             'unbound_from': lineEdit_unbound_from,
                             'unbound_to': lineEdit_unbound_to,
                             'radio_button': ''}
                item_now = place_dic[place_id]

                if place_id == 'range_from' or place_id == 'range_to':
                    global vis_adv_ranges_g

                    vis_adv_ranges_g = item_now.text()
                    validate_num(vis_adv_ranges, '')
                    if item_now.text() != vis_adv_ranges_g:
                        item_now.setText(vis_adv_ranges_g)

                elif place_id == 'unbound_from' or place_id == 'unbound_to':
                    if len(str(item_now.text())) >= 1:
                        try:
                            if type(int(item_now.text())) == int:
                                int_type = int(item_now.text())
                                if p_o_n == 0:
                                    if int_type == 0:
                                        item_now.setText('0')
                                    elif 0 < int_type <= 100:
                                        item_now.setText(str(int_type))
                                    elif int_type > 100:
                                        item_now.setText('100')
                                elif p_o_n == 1:
                                    if 0 < int_type:
                                        item_now.setText(str(int_type))
                        except:
                            item_now.setText('')

                elif place_id == 'radio_button':
                    if p_o_n == 0:
                        try:
                            if int(place_dic['unbound_from'].text()) > 100:
                                place_dic['unbound_from'].setText('100')
                        except:
                            pass
                        try:
                            if int(place_dic['unbound_to'].text()) > 100:
                                place_dic['unbound_to'].setText('100')
                        except:
                            pass

                update_add_button()

            def highlight_criteria(a):
                criteria_list_var = listWidget_list_criteria
                criteria_list = [criteria_list_var.item(index).text() for index in range(criteria_list_var.count())]
                selection = criteria_list[a]
                highlight_item = selection + ' - Highlighted'
                check_highlight = len([i.endswith('Highlighted') for i in criteria_list]) == 0
                highlight_cleaned = []

                for item in criteria_list:
                    if item.endswith('Highlighted'):
                        highlight_cleaned.append(item[:-14])
                    else:
                        highlight_cleaned.append(item)

                if not selection.endswith('Highlighted'):
                    if not check_highlight:
                        highlight_cleaned[a] = highlight_item
                        criteria_list_var.clear()
                        for h in highlight_cleaned:
                            criteria_list_var.insertItem(len(criteria_list), h)
                            pushButton_highlight_parameter.setDisabled(True)

                    elif check_highlight:
                        criteria_list[a] = highlight_item
                        criteria_list_var.clear()
                        for h in criteria_list:
                            criteria_list_var.insertItem(len(criteria_list), h)
                            pushButton_highlight_parameter.setDisabled(True)

                if selection.endswith('Highlighted'):
                    criteria_list_var.clear()
                    for h in highlight_cleaned:
                        criteria_list_var.insertItem(len(criteria_list), h)
                        pushButton_highlight_parameter.setText("Highlight")
                        pushButton_highlight_parameter.setDisabled(True)

                update_buttons('')

            def update_buttons(c_row):
                criteria_list_var = listWidget_list_criteria
                criteria_list = [criteria_list_var.item(index).text() for index in range(criteria_list_var.count())]
                selection = criteria_list[c_row] if c_row != '' else None

                if selection != None:
                    check_highlight = len([i.endswith('Highlighted') for i in criteria_list]) == 0
                    pushButton_delete_parameter.setEnabled(True)

                    if not check_highlight:
                        if not selection.endswith('Highlighted'):
                            pushButton_highlight_parameter.setEnabled(True)

                        elif selection.endswith('Highlighted'):
                            pushButton_highlight_parameter.setText("Un-highlight")
                            pushButton_highlight_parameter.setEnabled(True)

                    if not selection.endswith('Highlighted'):
                        pushButton_highlight_parameter.setText("Highlight")
                        pushButton_highlight_parameter.setEnabled(True)

                else:
                    pushButton_delete_parameter.setDisabled(True)
                    pushButton_highlight_parameter.setDisabled(True)

                if criteria_list != []:
                    buttonBox.setEnabled(True)
                    pushButton_reset_parameter.setEnabled(True)

                elif criteria_list == []:
                    buttonBox.setDisabled(True)
                    pushButton_reset_parameter.setDisabled(True)

            def add_criteria():
                item = get_add_criteria()[1]
                listWidget_list_criteria.insertItem(0, item)

                pushButton_add_parameter.setDisabled(True)
                update_buttons('')

            def set_ok():
                global advanced_criteria, advanced_criteria_text

                def convert_adv_criteria(d):
                    l_dummy = [None, None, None, None, None, False]

                    if 'Highlighted' in d:
                        l_dummy[5] = True

                    if '|' in d:
                        ran_un = d.split('|')
                        ran = ran_un[0][1:].split()
                        un = ran_un[1].split()
                        l_dummy[0] = int(ran[3]) if ran[3] != 'Min.' else None
                        l_dummy[1] = int(ran[5]) if ran[5] != 'Max.' else None
                        l_dummy[2] = int(un[3].strip('%')) if un[3] != 'Min.' else None
                        l_dummy[3] = int(un[5].strip('%')) if un[5] != 'Max.' else None
                        l_dummy[4] = un[1][:-1]

                    else:
                        ran = d[3:].split()
                        l_dummy[0] = int(ran[3]) if ran[3] != 'Min.' else None
                        l_dummy[1] = int(ran[5]) if ran[5] != 'Max.' else None

                    return l_dummy

                list_var = listWidget_list_criteria
                get_all = [list_var.item(index).text() for index in range(list_var.count())]

                del advanced_criteria_text[:]
                del advanced_criteria[:]

                adv_unsorted = []
                for c in get_all:
                    advanced_criteria_text.append(c)
                    cleaned = convert_adv_criteria(c)
                    adv_unsorted.append(cleaned)

                advanced_criteria = sorted(adv_unsorted, key=lambda x: x[5])
                self.checkBox_visual_advanced.setChecked(True)

            def delete_selected(c_row):
                criteria_list_var = listWidget_list_criteria
                criteria_list = [criteria_list_var.item(index).text() for index in range(criteria_list_var.count())]
                del criteria_list[c_row]
                listWidget_list_criteria.clear()
                for item in reversed(criteria_list):
                    listWidget_list_criteria.insertItem(0, item)

                update_add_button()
                update_buttons('')

            def reset_criteria_list():
                global advanced_criteria, advanced_criteria_text

                def get_user_consent_adv_run(consent_for):
                    consent_check = QMessageBox()
                    con_dic = {'reset_criteria_list': 'Resetting will remove all created criteria. '
                                                      'This action cannot be undone!' + '\n' + '\n'
                                                                                               'Click "Yes" to reset for a fresh start\n'
                                                                                               'Click "No" to return back', }

                    consent_check.setWindowTitle('Warning!')
                    consent_check.setText(con_dic[consent_for])
                    consent_check.setStandardButtons(
                        QMessageBox.Yes | QMessageBox.No)
                    consent = consent_check.exec_()

                    return consent

                consent = get_user_consent_adv_run('reset_criteria_list')

                if consent == 16384:
                    listWidget_list_criteria.clear()

                    advanced_criteria, advanced_criteria_text = [], []
                    self.checkBox_visual_advanced.setChecked(False)

                    update_buttons('')
                    update_add_button()

            def set_radio_status():
                global p_o_n
                checked_item = {1: radioButton_number.isChecked(),
                                0: radioButton_percentage.isChecked()}
                a = [i for i, status in checked_item.items() if status]
                p_o_n = a[0]
                validate_filter_input('radio_button')
                update_add_button()

            def set_existing_parameters():
                if advanced_criteria_text != []:
                    for item in reversed(advanced_criteria_text):
                        listWidget_list_criteria.insertItem(0, item)

                    update_buttons('')

            visual_adv_dialog = QDialog()
            visual_adv_dialog.setObjectName("Dialog")
            visual_adv_dialog.setEnabled(True)
            visual_adv_dialog.resize(800, 350)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(visual_adv_dialog.sizePolicy().hasHeightForWidth())
            visual_adv_dialog.setSizePolicy(sizePolicy)
            visual_adv_dialog.setMaximumSize(QtCore.QSize(900, 405))
            gridLayout_3 = QtWidgets.QGridLayout(visual_adv_dialog)
            gridLayout_3.setObjectName("gridLayout_3")
            frame_5 = QtWidgets.QFrame(visual_adv_dialog)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_5.sizePolicy().hasHeightForWidth())
            frame_5.setSizePolicy(sizePolicy)
            frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_5.setObjectName("frame_5")
            gridLayout = QtWidgets.QGridLayout(frame_5)
            gridLayout.setContentsMargins(20, 20, 20, 20)
            gridLayout.setObjectName("gridLayout")
            gridLayout_2 = QtWidgets.QGridLayout()
            gridLayout_2.setObjectName("gridLayout_2")
            horizontalLayout = QtWidgets.QHBoxLayout()
            horizontalLayout.setObjectName("horizontalLayout")
            label_ssDNA_range = QtWidgets.QLabel(frame_5)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(label_ssDNA_range.sizePolicy().hasHeightForWidth())
            label_ssDNA_range.setSizePolicy(sizePolicy)
            label_ssDNA_range.setObjectName("label_ssDNA_range")
            horizontalLayout.addWidget(label_ssDNA_range)
            lineEdit_range_from = QtWidgets.QLineEdit(frame_5)
            font = QtGui.QFont()
            font.setPointSize(9)
            lineEdit_range_from.setFont(font)
            lineEdit_range_from.setObjectName("lineEdit_range_from")
            horizontalLayout.addWidget(lineEdit_range_from)
            label_2 = QtWidgets.QLabel(frame_5)
            label_2.setObjectName("label_2")
            horizontalLayout.addWidget(label_2)
            lineEdit_range_to = QtWidgets.QLineEdit(frame_5)
            font = QtGui.QFont()
            font.setPointSize(9)
            lineEdit_range_to.setFont(font)
            lineEdit_range_to.setObjectName("lineEdit_range_to")
            horizontalLayout.addWidget(lineEdit_range_to)
            label_5 = QtWidgets.QLabel(frame_5)
            font = QtGui.QFont()
            font.setPointSize(12)
            font.setBold(True)
            font.setUnderline(False)
            font.setWeight(75)
            font.setKerning(True)
            label_5.setFont(font)
            label_5.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
            label_5.setObjectName("label_5")
            horizontalLayout.addWidget(label_5)
            label_unbound = QtWidgets.QLabel(frame_5)
            label_unbound.setObjectName("label_unbound")
            horizontalLayout.addWidget(label_unbound)
            radioButton_percentage = QtWidgets.QRadioButton(frame_5)
            radioButton_percentage.setObjectName("radioButton_percentage")
            horizontalLayout.addWidget(radioButton_percentage)
            radioButton_number = QtWidgets.QRadioButton(frame_5)
            radioButton_number.setObjectName("radioButton_number")
            horizontalLayout.addWidget(radioButton_number)
            lineEdit_unbound_from = QtWidgets.QLineEdit(frame_5)
            font = QtGui.QFont()
            font.setPointSize(9)
            lineEdit_unbound_from.setFont(font)
            lineEdit_unbound_from.setObjectName("lineEdit_unbound_from")
            horizontalLayout.addWidget(lineEdit_unbound_from)
            label_3 = QtWidgets.QLabel(frame_5)
            label_3.setObjectName("label_3")
            horizontalLayout.addWidget(label_3)
            lineEdit_unbound_to = QtWidgets.QLineEdit(frame_5)
            font = QtGui.QFont()
            font.setPointSize(9)
            lineEdit_unbound_to.setFont(font)
            lineEdit_unbound_to.setObjectName("lineEdit_unbound_to")
            horizontalLayout.addWidget(lineEdit_unbound_to)
            pushButton_add_parameter = QtWidgets.QPushButton(frame_5)
            pushButton_add_parameter.setObjectName("pushButton_add_parameter")
            horizontalLayout.addWidget(pushButton_add_parameter)
            gridLayout_2.addLayout(horizontalLayout, 0, 0, 1, 1)
            frame = QtWidgets.QFrame(frame_5)
            frame.setMinimumSize(QtCore.QSize(0, 15))
            frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame.setFrameShadow(QtWidgets.QFrame.Raised)
            frame.setObjectName("frame")
            gridLayout_2.addWidget(frame, 1, 0, 1, 1)
            gridLayout_parameters = QtWidgets.QGridLayout()
            gridLayout_parameters.setContentsMargins(-1, 0, -1, -1)
            gridLayout_parameters.setObjectName("gridLayout_parameters")
            listWidget_list_criteria = QtWidgets.QListWidget(frame_5)
            listWidget_list_criteria.setMinimumSize(QtCore.QSize(0, 150))
            listWidget_list_criteria.setObjectName("listWidget_list_criteria")
            gridLayout_parameters.addWidget(listWidget_list_criteria, 0, 0, 7, 1)
            pushButton_highlight_parameter = QtWidgets.QPushButton(frame_5)
            pushButton_highlight_parameter.setObjectName("pushButton_highlight_parameter")
            gridLayout_parameters.addWidget(pushButton_highlight_parameter, 1, 1, 1, 1)
            frame_6 = QtWidgets.QFrame(frame_5)
            frame_6.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_6.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_6.setObjectName("frame_6")
            gridLayout_parameters.addWidget(frame_6, 6, 1, 1, 1)
            pushButton_reset_parameter = QtWidgets.QPushButton(frame_5)
            pushButton_reset_parameter.setObjectName("pushButton_reset_parameter")
            gridLayout_parameters.addWidget(pushButton_reset_parameter, 5, 1, 1, 1)
            frame_4 = QtWidgets.QFrame(frame_5)
            frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_4.setObjectName("frame_4")
            gridLayout_parameters.addWidget(frame_4, 4, 1, 1, 1)
            frame_3 = QtWidgets.QFrame(frame_5)
            frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_3.setObjectName("frame_3")
            gridLayout_parameters.addWidget(frame_3, 2, 1, 1, 1)
            pushButton_delete_parameter = QtWidgets.QPushButton(frame_5)
            pushButton_delete_parameter.setObjectName("pushButton_delete_parameter")
            gridLayout_parameters.addWidget(pushButton_delete_parameter, 3, 1, 1, 1)
            frame_7 = QtWidgets.QFrame(frame_5)
            frame_7.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_7.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_7.setObjectName("frame_7")
            gridLayout_parameters.addWidget(frame_7, 0, 1, 1, 1)
            gridLayout_2.addLayout(gridLayout_parameters, 2, 0, 1, 1)
            frame_2 = QtWidgets.QFrame(frame_5)
            frame_2.setMinimumSize(QtCore.QSize(0, 15))
            frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_2.setObjectName("frame_2")
            gridLayout_2.addWidget(frame_2, 3, 0, 1, 1)
            horizontalLayout_ok_cancel = QtWidgets.QHBoxLayout()
            horizontalLayout_ok_cancel.setObjectName("horizontalLayout_ok_cancel")
            buttonBox = QtWidgets.QDialogButtonBox(frame_5)
            buttonBox.setOrientation(QtCore.Qt.Horizontal)
            buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Ok)
            buttonBox.setObjectName("buttonBox")
            horizontalLayout_ok_cancel.addWidget(buttonBox)
            pushButton_visual_cancel = QtWidgets.QPushButton(frame_5)
            pushButton_visual_cancel.setObjectName("pushButton_visual_cancel")
            horizontalLayout_ok_cancel.addWidget(pushButton_visual_cancel)
            gridLayout_2.addLayout(horizontalLayout_ok_cancel, 4, 0, 1, 1)
            gridLayout.addLayout(gridLayout_2, 0, 0, 1, 1)
            gridLayout_3.addWidget(frame_5, 0, 0, 1, 1)

            def retranslateUi(Dialog):
                _translate = QtCore.QCoreApplication.translate
                Dialog.setWindowTitle(_translate("Dialog", "Advanced filtering"))
                label_ssDNA_range.setText(_translate("Dialog", "ssDNA(s) range from"))
                label_2.setText(_translate("Dialog", "to"))
                label_5.setText(_translate("Dialog", "|"))
                label_unbound.setText(_translate("Dialog", "Unbound"))
                radioButton_percentage.setText(_translate("Dialog", "%"))
                radioButton_number.setText(_translate("Dialog", "n"))
                label_3.setText(_translate("Dialog", "to"))
                pushButton_add_parameter.setText(_translate("Dialog", "Add"))
                pushButton_highlight_parameter.setText(_translate("Dialog", "Highlight"))
                pushButton_reset_parameter.setText(_translate("Dialog", "Reset all"))
                pushButton_delete_parameter.setText(_translate("Dialog", "Delete"))
                pushButton_visual_cancel.setText(_translate("Dialog", "Cancel"))

            retranslateUi(visual_adv_dialog)
            buttonBox.accepted.connect(visual_adv_dialog.accept)
            buttonBox.rejected.connect(visual_adv_dialog.reject)

            def close_dialog():
                visual_adv_dialog.reject()

            pushButton_add_parameter.clicked.connect(add_criteria)
            radioButton_number.clicked.connect(set_radio_status)
            radioButton_percentage.clicked.connect(set_radio_status)
            radioButton_percentage.setChecked(True)
            pushButton_add_parameter.setDisabled(True)
            lineEdit_range_from.textChanged.connect(lambda: validate_filter_input('range_from'))
            lineEdit_range_to.textChanged.connect(lambda: validate_filter_input('range_to'))
            lineEdit_unbound_from.textChanged.connect(lambda: validate_filter_input('unbound_from'))
            lineEdit_unbound_to.textChanged.connect(lambda: validate_filter_input('unbound_to'))
            listWidget_list_criteria.itemClicked.connect(lambda: update_buttons(listWidget_list_criteria.currentRow()))
            pushButton_highlight_parameter.clicked.connect(lambda:
                                                           highlight_criteria(listWidget_list_criteria.currentRow()))
            buttonBox.setDisabled(True)
            pushButton_delete_parameter.setDisabled(True)
            pushButton_highlight_parameter.setDisabled(True)
            pushButton_reset_parameter.setDisabled(True)
            pushButton_delete_parameter.clicked.connect(lambda: delete_selected(listWidget_list_criteria.currentRow()))
            pushButton_reset_parameter.clicked.connect(reset_criteria_list)
            pushButton_visual_cancel.clicked.connect(close_dialog)
            set_existing_parameters() if advanced_criteria_text != [] else None
            visualize_advanced = visual_adv_dialog.exec_()

            if visualize_advanced == 1:
                set_ok()

        def edit_seq(c_row):
            to_be_edited_seq = ssdna_comp_list[c_row][2]
            to_be_edited_amount = str(ssdna_comp_list[c_row][1])
            to_be_edited_type = ssdna_comp_list[c_row][0]

            def update_edit_buttons():
                seqs_in = lineEdit_write_edit_sequence.text()
                amount_in = lineEdit_write_edit_amount.text()

                if to_be_edited_type == 1:
                    exist_items = [e[2] for e in ssdna_comp_list]
                    del exist_items[c_row]

                    if seqs_in != '' and amount_in != '':
                        if seqs_in != to_be_edited_seq or amount_in != to_be_edited_amount:
                            if len(seqs_in) >= 2 and len(amount_in) >= 1:
                                if seqs_in not in exist_items:
                                    buttonBox_write_edit.setEnabled(True)
                                else:
                                    buttonBox_write_edit.setDisabled(True)
                            else:
                                buttonBox_write_edit.setDisabled(True)
                        else:
                            buttonBox_write_edit.setDisabled(True)
                    else:
                        buttonBox_write_edit.setDisabled(True)

                elif to_be_edited_type == 2:
                    if amount_in != '':
                        if amount_in != to_be_edited_amount:
                            if len(amount_in) >= 1:
                                buttonBox_write_edit.setEnabled(True)
                            else:
                                buttonBox_write_edit.setDisabled(True)
                        else:
                            buttonBox_write_edit.setDisabled(True)
                    else:
                        buttonBox_write_edit.setDisabled(True)
                else:
                    buttonBox_write_edit.setDisabled(True)

            def save_edited_seq(item_row, edited_seq, edited_seq_amount):
                global ssdna_comp_list

                seq_syn = create_bngl_sequence(edited_seq) if to_be_edited_seq != edited_seq \
                    else ssdna_comp_list[item_row][3]
                get_edited = [to_be_edited_type, int(edited_seq_amount), edited_seq, seq_syn]
                ssdna_comp_list[item_row] = get_edited

                update_created_list()

                return True

            def set_validated(call_from):
                global write_edited_amount_g, write_edited_seq_g
                if call_from == write_edit_amount:
                    write_edited_amount_g = lineEdit_write_edit_amount.text()
                    validate_num(call_from, '')
                    if lineEdit_write_edit_amount.text() != write_edited_amount_g:
                        lineEdit_write_edit_amount.setText(write_edited_amount_g)

                elif call_from == write_edit_seq:
                    write_edited_seq_g = lineEdit_write_edit_sequence.text()
                    validate_char(call_from)
                    if lineEdit_write_edit_sequence.text() != write_edited_seq_g:
                        lineEdit_write_edit_sequence.setText(write_edited_seq_g)

                update_edit_buttons()

            write_edit_seq_dialog = QDialog()
            write_edit_seq_dialog.setObjectName("Dialog")
            write_edit_seq_dialog.resize(926, 162)
            write_edit_seq_dialog.setMaximumSize(QtCore.QSize(1157, 202))
            gridLayout_2 = QtWidgets.QGridLayout(write_edit_seq_dialog)
            gridLayout_2.setObjectName("gridLayout_2")
            frame_2 = QtWidgets.QFrame(write_edit_seq_dialog)
            frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_2.setObjectName("frame_2")
            gridLayout = QtWidgets.QGridLayout(frame_2)
            gridLayout.setObjectName("gridLayout")
            verticalLayout = QtWidgets.QVBoxLayout()
            verticalLayout.setObjectName("verticalLayout")
            gridLayout_write_edit_sequence_entry = QtWidgets.QGridLayout()
            gridLayout_write_edit_sequence_entry.setObjectName("gridLayout_write_edit_sequence_entry")
            frame_32 = QtWidgets.QFrame(frame_2)
            frame_32.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_32.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_32.setObjectName("frame_32")
            gridLayout_write_edit_sequence_entry.addWidget(frame_32, 2, 0, 1, 1)
            frame_8 = QtWidgets.QFrame(frame_2)
            frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_8.setObjectName("frame_8")
            gridLayout_write_edit_sequence_entry.addWidget(frame_8, 0, 2, 1, 1)
            frame_30 = QtWidgets.QFrame(frame_2)
            frame_30.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_30.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_30.setObjectName("frame_30")
            gridLayout_write_edit_sequence_entry.addWidget(frame_30, 2, 2, 1, 1)
            frame_7 = QtWidgets.QFrame(frame_2)
            frame_7.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_7.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_7.setObjectName("frame_7")
            gridLayout_write_edit_sequence_entry.addWidget(frame_7, 0, 1, 1, 1)
            label_write_edit_amount = QtWidgets.QLabel(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(label_write_edit_amount.sizePolicy().hasHeightForWidth())
            label_write_edit_amount.setSizePolicy(sizePolicy)
            label_write_edit_amount.setAlignment(QtCore.Qt.AlignCenter)
            label_write_edit_amount.setObjectName("label_write_edit_amount")
            gridLayout_write_edit_sequence_entry.addWidget(label_write_edit_amount, 0, 4, 1, 1)
            label_write_edit_allowed_num = QtWidgets.QLabel(frame_2)
            font = QtGui.QFont()
            font.setPointSize(8)
            font.setItalic(True)
            label_write_edit_allowed_num.setFont(font)
            label_write_edit_allowed_num.setAlignment(QtCore.Qt.AlignCenter)
            label_write_edit_allowed_num.setObjectName("label_write_edit_allowed_num")
            gridLayout_write_edit_sequence_entry.addWidget(label_write_edit_allowed_num, 2, 4, 1, 1)
            label_write_edit_sequence = QtWidgets.QLabel(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(label_write_edit_sequence.sizePolicy().hasHeightForWidth())
            label_write_edit_sequence.setSizePolicy(sizePolicy)
            label_write_edit_sequence.setAlignment(
                QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            label_write_edit_sequence.setObjectName("label_write_edit_sequence")
            gridLayout_write_edit_sequence_entry.addWidget(label_write_edit_sequence, 0, 0, 1, 1)
            label_write_edit_allowed_char = QtWidgets.QLabel(frame_2)
            font = QtGui.QFont()
            font.setPointSize(8)
            font.setItalic(True)
            label_write_edit_allowed_char.setFont(font)
            label_write_edit_allowed_char.setAlignment(QtCore.Qt.AlignCenter)
            label_write_edit_allowed_char.setObjectName("label_write_edit_allowed_char")
            gridLayout_write_edit_sequence_entry.addWidget(label_write_edit_allowed_char, 2, 1, 1, 1)
            lineEdit_write_edit_amount = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_edit_amount.sizePolicy().hasHeightForWidth())
            lineEdit_write_edit_amount.setSizePolicy(sizePolicy)
            font = QtGui.QFont()
            font.setPointSize(9)
            lineEdit_write_edit_amount.setFont(font)
            lineEdit_write_edit_amount.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_edit_amount.setObjectName("lineEdit_write_edit_amount")
            gridLayout_write_edit_sequence_entry.addWidget(lineEdit_write_edit_amount, 1, 4, 1, 1)
            frame_3 = QtWidgets.QFrame(frame_2)
            frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_3.setObjectName("frame_3")
            gridLayout_write_edit_sequence_entry.addWidget(frame_3, 0, 3, 1, 1)
            lineEdit_write_edit_sequence = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_edit_sequence.sizePolicy().hasHeightForWidth())
            lineEdit_write_edit_sequence.setSizePolicy(sizePolicy)
            font = QtGui.QFont()
            font.setPointSize(9)
            lineEdit_write_edit_sequence.setFont(font)
            lineEdit_write_edit_sequence.setObjectName("lineEdit_write_edit_sequence")
            gridLayout_write_edit_sequence_entry.addWidget(lineEdit_write_edit_sequence, 1, 0, 1, 4)
            frame_4 = QtWidgets.QFrame(frame_2)
            frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_4.setObjectName("frame_4")
            gridLayout_write_edit_sequence_entry.addWidget(frame_4, 2, 3, 1, 1)
            verticalLayout.addLayout(gridLayout_write_edit_sequence_entry)
            frame = QtWidgets.QFrame(frame_2)
            frame.setMinimumSize(QtCore.QSize(0, 15))
            frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame.setFrameShadow(QtWidgets.QFrame.Raised)
            frame.setObjectName("frame")
            verticalLayout.addWidget(frame)
            horizontalLayout_save_cancel = QtWidgets.QHBoxLayout()
            horizontalLayout_save_cancel.setObjectName("horizontalLayout_save_cancel")
            buttonBox_write_edit = QtWidgets.QDialogButtonBox(frame_2)
            buttonBox_write_edit.setOrientation(QtCore.Qt.Horizontal)
            buttonBox_write_edit.setStandardButtons(QtWidgets.QDialogButtonBox.Save)
            buttonBox_write_edit.setObjectName("buttonBox_write_edit")
            horizontalLayout_save_cancel.addWidget(buttonBox_write_edit)
            pushButton_cancel = QtWidgets.QPushButton(frame_2)
            pushButton_cancel.setObjectName("pushButton_cancel")
            horizontalLayout_save_cancel.addWidget(pushButton_cancel)
            verticalLayout.addLayout(horizontalLayout_save_cancel)
            gridLayout.addLayout(verticalLayout, 0, 0, 1, 1)
            gridLayout_2.addWidget(frame_2, 1, 0, 1, 1)
            buttonBox_write_edit.setDisabled(True)

            def retranslateUi(Dialog):
                _translate = QtCore.QCoreApplication.translate
                Dialog.setWindowTitle(_translate("Dialog", "Edit sequence"))
                label_write_edit_amount.setText(_translate("Dialog", "Amount"))
                label_write_edit_allowed_num.setText(_translate("Dialog", "Numericals only"))
                label_write_edit_sequence.setText(_translate("Dialog", "Edit ssDNA sequence"))
                label_write_edit_allowed_char.setText(_translate("Dialog", "Allowed characters (A, T, C, G)"))
                pushButton_cancel.setText(_translate("Dialog", "Cancel"))

            def close_dialog():
                write_edit_seq_dialog.reject()

            retranslateUi(write_edit_seq_dialog)
            buttonBox_write_edit.accepted.connect(write_edit_seq_dialog.accept)
            buttonBox_write_edit.rejected.connect(write_edit_seq_dialog.reject)
            write_edited_amount_g = lineEdit_write_edit_amount

            lineEdit_write_edit_sequence.setText(to_be_edited_seq)
            lineEdit_write_edit_amount.setText(to_be_edited_amount)
            lineEdit_write_edit_amount.textChanged.connect(lambda: set_validated(write_edit_amount))
            lineEdit_write_edit_sequence.textChanged.connect(lambda: set_validated(write_edit_seq))
            pushButton_cancel.clicked.connect(close_dialog)
            lineEdit_write_edit_sequence.setDisabled(True) if ssdna_comp_list[c_row][0] == 2 \
                else lineEdit_write_edit_sequence.setEnabled(True)
            run_write_edit = write_edit_seq_dialog.exec_()

            if run_write_edit == 1:
                edited_seq = lineEdit_write_edit_sequence.text()
                edited_amount = lineEdit_write_edit_amount.text()
                save_seq = save_edited_seq(c_row, edited_seq, edited_amount)

                if save_seq:
                    self.label_messages.setText('Sequence edited to successfully ✔')
                    validate_and_update_buttons(write_tab)

        def get_k_values(f_p):
            data, q = {}, []
            with open(f_p, 'r') as f:
                begin_k, end_k = 'begin parameters', 'end parameters'
                l = f.readlines()
                for i in l:
                    check_syn = i.rstrip('\n').strip()

                    if begin_k not in q:
                        if check_syn != begin_k:
                            continue

                        elif check_syn == begin_k:
                            q.append(begin_k)

                    elif begin_k in q:
                        if check_syn.startswith('k'):
                            check_k = check_syn.split(' ')

                            if len(check_k) == 2:
                                if check_k[0][0] == 'k' and check_k[0][1:].isdigit():
                                    if check_k[1].isdigit():
                                        k_v = 'k' + str(check_k[0][1:])
                                        data.update({k_v: check_k[-1]})

                                    else:
                                        try:
                                            try_float = float(check_k[1])
                                            k_v = 'k' + str(check_k[0][1:])
                                            data.update({k_v: str(try_float)})
                                        except:
                                            pass
                            else:
                                continue

                        elif check_syn.startswith('Temp'):
                            check_k = check_syn.split(' ')
                            if check_k[0] == 'Temp' and check_k[-1].isdigit():
                                data.update({check_k[0]: check_k[-1]})

                            else:
                                continue

                        elif check_syn == end_k:
                            break
                        else:
                            continue
            return data

        def run_activate_run_advanced():
            if self.checkBox_run_advanced.isChecked():
                check_is_adv_data_1 = sum([1 for i in run_temper_data_g if i != '']) == 4
                check_is_adv_data_2 = sum([1 for i in run_time_data_g if i != '']) == 4

                if run_kt_list_g != {}:
                    if run_annealing_g and check_is_adv_data_1 and check_is_adv_data_2:
                        self.label_messages.setText('')
                        self.lineEdit_run_start_time.setDisabled(True)
                        self.lineEdit_run_n_dumps.setDisabled(True)
                        self.lineEdit_run_sim_end.setDisabled(True)

                    else:
                        self.lineEdit_run_start_time.setEnabled(True)
                        self.lineEdit_run_n_dumps.setEnabled(True)
                        self.lineEdit_run_sim_end.setEnabled(True)

                else:
                    self.checkBox_run_advanced.setChecked(False)
                    self.lineEdit_run_start_time.setEnabled(True)
                    self.lineEdit_run_n_dumps.setEnabled(True)
                    self.lineEdit_run_sim_end.setEnabled(True)
                    self.label_messages.setText('Advanced parameters not entered in advanced menu ✘')

            elif not self.checkBox_run_advanced.isChecked():
                self.lineEdit_run_start_time.setEnabled(True)
                self.lineEdit_run_n_dumps.setEnabled(True)
                self.lineEdit_run_sim_end.setEnabled(True)
                self.label_messages.setText('')

        def run_advanced():
            source_file = self.lineEdit_run_browse_source.text()

            auto_cal_list = ['run_start_time', 'run_time_per_stage', 'run_dumps_per_stage',
                             'run_start_temp', 'run_d_temp', 'run_end_temp']

            def set_previous():
                lineEdit_run_start_time.setText(run_time_data_g[0])
                lineEdit_run_time_per_stage.setText(run_time_data_g[1])
                lineEdit_run_dumps_per_stage.setText(run_time_data_g[2])

                lineEdit_run_start_temp.setText(run_temper_data_g[0])
                lineEdit_run_d_temp.setText(run_temper_data_g[1])
                lineEdit_run_end_temp.setText(run_temper_data_g[2])

            def update_save_button():

                if checkBox_annealing.isChecked():
                    all_cells_check = sum([1 for i in kt_n_ann_cells.values() if i.text().endswith('.')]) == 0
                    annealing_data = [lineEdit_run_start_temp.text(), lineEdit_run_d_temp.text(),
                                      lineEdit_run_end_temp.text(),
                                      lineEdit_run_start_time.text(), lineEdit_run_time_per_stage.text(),
                                      lineEdit_run_dumps_per_stage.text()]

                    data_exists = sum([1 for i in annealing_data if i.strip() != '']) == 6

                    if data_exists and all_cells_check:
                        buttonBox_write_edit.setEnabled(True)

                    elif not data_exists or not all_cells_check:
                        buttonBox_write_edit.setDisabled(True)

                if not checkBox_annealing.isChecked():
                    kt_cells_check = sum([1 for i in kt_n_ann_cells.items() if i[0] not in auto_cal_list
                                                                            and i[1].text().endswith('.')]) == 0
                    if kt_cells_check:
                        buttonBox_write_edit.setEnabled(True)

                    elif not kt_cells_check:
                        buttonBox_write_edit.setDisabled(True)

            def set_k_values(f_path, call_from):

                kt_cells = {'Temp': [lineEdit_run_temp_p, lineEdit_run_temp_vdna],
                            'k1': [lineEdit_run_p_k1, lineEdit_run_vdna_k1],
                            'k2': [lineEdit_run_p_k2, lineEdit_run_vdna_k2],
                            'k3': [lineEdit_run_p_k3, lineEdit_run_vdna_k3],
                            'k4': [lineEdit_run_p_k4, lineEdit_run_vdna_k4],
                            'k5': [lineEdit_run_p_k5, lineEdit_run_vdna_k5],
                            'k6': [lineEdit_run_p_k6, lineEdit_run_vdna_k6],
                            'k7': [lineEdit_run_p_k7, lineEdit_run_vdna_k7],
                            'k8': [lineEdit_run_p_k8, lineEdit_run_vdna_k8],
                            'k9': [lineEdit_run_p_k9, lineEdit_run_vdna_k9],
                            'k10': [lineEdit_run_p_k10, lineEdit_run_vdna_k10],
                            'k11': [lineEdit_run_p_k11, lineEdit_run_vdna_k11]}

                def set_file_run_parameters(a):
                    k_data_r = get_k_values('system_files/reference_files/ref_file.bngl')

                    if a == 'get':
                        k_data_p = get_k_values(f_path) if run_kt_list_g == {} else run_kt_list_g
                        for ii in kt_cells.items():
                            if ii[0] in k_data_p:
                                ii[1][0].setText(k_data_p[ii[0]])
                            else:
                                ii[1][0].setText('')

                            if ii[0] in k_data_r:
                                ii[1][1].setText(k_data_r[ii[0]])
                            else:
                                ii[1][1].setText('')

                    elif a == 'set_d':
                        for ii in kt_cells.items():
                            if ii[0] in k_data_r:
                                ii[1][0].setText(k_data_r[ii[0]])
                            else:
                                ii[1][0].setText('')

                    elif a == 'set_o':
                        k_data_p = get_k_values(f_path)
                        for ii in kt_cells.items():
                            if ii[0] in k_data_p:
                                ii[1][0].setText(k_data_p[ii[0]])
                            else:
                                ii[1][0].setText('')

                set_file_run_parameters(call_from)

            def save_rnf_data():
                global run_kt_list_g, run_temper_data_g, run_time_data_g

                k_list_values = {'Temp': lineEdit_run_temp_p.text(),
                                 'k1': lineEdit_run_p_k1.text(), 'k2': lineEdit_run_p_k2.text(),
                                 'k3': lineEdit_run_p_k3.text(), 'k4': lineEdit_run_p_k4.text(),
                                 'k5': lineEdit_run_p_k5.text(), 'k6': lineEdit_run_p_k6.text(),
                                 'k7': lineEdit_run_p_k7.text(), 'k8': lineEdit_run_p_k8.text(),
                                 'k9': lineEdit_run_p_k9.text(), 'k10': lineEdit_run_p_k10.text(),
                                 'k11': lineEdit_run_p_k11.text()}

                del_k = []
                for key, val in k_list_values.items():
                    if val == '':
                        del_k.append(key)

                for k in del_k:
                    del k_list_values[k]

                run_kt_list_g = k_list_values

                run_time_data_g = [lineEdit_run_start_time.text(), lineEdit_run_time_per_stage.text(),
                                   lineEdit_run_dumps_per_stage.text(), lineEdit_run_total_time.text()]

                run_temper_data_g = [lineEdit_run_start_temp.text(), lineEdit_run_d_temp.text(),
                                     lineEdit_run_end_temp.text(), lineEdit_run_n_stages.text()]

            def auto_cal():

                temp_data = [lineEdit_run_start_temp.text(), lineEdit_run_d_temp.text(),
                             lineEdit_run_end_temp.text()]
                time_data = [lineEdit_run_start_time.text(),
                             lineEdit_run_time_per_stage.text(), lineEdit_run_total_time.text()]

                float_temp_list = []
                float_time_list = []
                try:
                    float_temp_list_1 = [float(i) if i != '' else None for i in temp_data]
                    float_time_list_1 = [float(i) if i != '' else None for i in time_data]
                    float_temp_list = float_temp_list_1
                    float_time_list = float_time_list_1

                except:
                    pass

                try:
                    reverse_if_temp = sorted([float_temp_list[0], float_temp_list[2]])[::-1]
                    reverse_if_temp.insert(1, float_temp_list[1])

                    float_temp_list = reverse_if_temp
                except:
                    pass

                try:
                    n_stages = ceil(((float_temp_list[0] - float_temp_list[2]) / float_temp_list[1]) + 1)
                    lineEdit_run_n_stages.setText(str(n_stages))
                except:
                    lineEdit_run_n_stages.setText('')

                try:
                    total_time = round(float_time_list[1] * int(lineEdit_run_n_stages.text()), 4)
                    lineEdit_run_total_time.setText(str(total_time))
                except:
                    lineEdit_run_total_time.setText('')

                update_save_button()

            def validate_cell_values(call_from):
                global write_n_run_kt_validated_g

                current_item = kt_n_ann_cells[call_from].text()

                write_n_run_kt_validated_g = current_item
                validate_num(write_n_run_kt_validate, '')

                if current_item != write_n_run_kt_validated_g:
                    kt_n_ann_cells[call_from].setText(write_n_run_kt_validated_g)
                    write_n_run_kt_validated_g = ''

                if call_from in auto_cal_list:
                    auto_cal()
                else:
                    update_save_button()

            def activate_annealing(call_from):
                annealing_data = [lineEdit_run_start_time, lineEdit_run_time_per_stage,
                                  lineEdit_run_dumps_per_stage, lineEdit_run_total_time,
                                  lineEdit_run_start_temp, lineEdit_run_d_temp,
                                  lineEdit_run_end_temp, lineEdit_run_n_stages]

                check_is_adv_data_1 = sum([1 for i in run_temper_data_g if i != '']) == 4
                check_is_adv_data_2 = sum([1 for i in run_time_data_g if i != '']) == 4

                if not call_from == 'u_clicked':
                    if run_annealing_g:
                        checkBox_annealing.setChecked(True)
                        for cell in annealing_data:
                            cell.setEnabled(True)
                        set_previous()

                    elif not run_annealing_g:
                        checkBox_annealing.setChecked(False)
                        for cell in annealing_data:
                            cell.setDisabled(True)
                        set_previous()

                if call_from == 'u_clicked':
                    if not checkBox_annealing.isChecked():
                        for cell in annealing_data:
                            cell.setDisabled(True)

                    elif checkBox_annealing.isChecked():
                        for cell in annealing_data:
                            cell.setEnabled(True)
                            set_previous()
                        if check_is_adv_data_1 and check_is_adv_data_2:
                            set_previous()

                update_save_button()

            run_advanced_dialog = QDialog()
            run_advanced_dialog.setObjectName("Dialog")
            run_advanced_dialog.resize(600, 770)
            run_advanced_dialog.setMaximumSize(QtCore.QSize(700, 810))
            gridLayout_2 = QtWidgets.QGridLayout(run_advanced_dialog)
            gridLayout_2.setObjectName("gridLayout_2")
            frame_2 = QtWidgets.QFrame(run_advanced_dialog)
            frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_2.setObjectName("frame_2")
            gridLayout = QtWidgets.QGridLayout(frame_2)
            gridLayout.setContentsMargins(20, 20, 20, 20)
            gridLayout.setObjectName("gridLayout")
            verticalLayout = QtWidgets.QVBoxLayout()
            verticalLayout.setObjectName("verticalLayout")
            gridLayout_run_parameters_entry = QtWidgets.QGridLayout()
            gridLayout_run_parameters_entry.setObjectName("gridLayout_run_parameters_entry")
            lineEdit_run_n_stages = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_n_stages.sizePolicy().hasHeightForWidth())
            lineEdit_run_n_stages.setSizePolicy(sizePolicy)
            lineEdit_run_n_stages.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_n_stages.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_n_stages.setReadOnly(True)
            lineEdit_run_n_stages.setObjectName("lineEdit_run_n_stages")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_n_stages, 18, 6, 1, 1)
            lineEdit_run_start_time = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_start_time.sizePolicy().hasHeightForWidth())
            lineEdit_run_start_time.setSizePolicy(sizePolicy)
            lineEdit_run_start_time.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_start_time.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_start_time.setObjectName("lineEdit_run_start_time")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_start_time, 21, 0, 2, 1)
            lineEdit_run_p_k9 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k9.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k9.setSizePolicy(sizePolicy)
            lineEdit_run_p_k9.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k9.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k9.setObjectName("lineEdit_run_p_k9")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k9, 9, 2, 1, 1)
            lineEdit_run_vdna_k5 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k5.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k5.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k5.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k5.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k5.setReadOnly(True)
            lineEdit_run_vdna_k5.setObjectName("lineEdit_run_vdna_k5")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k5, 5, 4, 1, 1)
            label_run_vdna_defaults = QtWidgets.QLabel(frame_2)
            label_run_vdna_defaults.setAlignment(QtCore.Qt.AlignCenter)
            label_run_vdna_defaults.setObjectName("label_run_vdna_defaults")
            gridLayout_run_parameters_entry.addWidget(label_run_vdna_defaults, 0, 4, 1, 1)
            lineEdit_run_p_k4 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k4.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k4.setSizePolicy(sizePolicy)
            lineEdit_run_p_k4.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k4.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k4.setObjectName("lineEdit_run_p_k4")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k4, 4, 2, 1, 1)
            lineEdit_run_vdna_k1 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k1.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k1.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k1.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k1.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k1.setReadOnly(True)
            lineEdit_run_vdna_k1.setObjectName("lineEdit_run_vdna_k1")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k1, 1, 4, 1, 1)
            lineEdit_run_p_k10 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k10.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k10.setSizePolicy(sizePolicy)
            lineEdit_run_p_k10.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k10.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k10.setObjectName("lineEdit_run_p_k10")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k10, 10, 2, 1, 1)
            label_run_k4 = QtWidgets.QLabel(frame_2)
            label_run_k4.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k4.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k4.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k4.setObjectName("label_run_k4")
            gridLayout_run_parameters_entry.addWidget(label_run_k4, 4, 0, 1, 1)
            lineEdit_run_p_k1 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k1.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k1.setSizePolicy(sizePolicy)
            lineEdit_run_p_k1.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k1.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k1.setObjectName("lineEdit_run_p_k1")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k1, 1, 2, 1, 1)
            label_run_k3 = QtWidgets.QLabel(frame_2)
            label_run_k3.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k3.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k3.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k3.setObjectName("label_run_k3")
            gridLayout_run_parameters_entry.addWidget(label_run_k3, 3, 0, 1, 1)
            lineEdit_run_vdna_k8 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k8.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k8.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k8.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k8.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k8.setReadOnly(True)
            lineEdit_run_vdna_k8.setObjectName("lineEdit_run_vdna_k8")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k8, 8, 4, 1, 1)
            lineEdit_run_p_k2 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k2.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k2.setSizePolicy(sizePolicy)
            lineEdit_run_p_k2.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k2.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k2.setObjectName("lineEdit_run_p_k2")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k2, 2, 2, 1, 1)
            lineEdit_run_p_k7 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k7.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k7.setSizePolicy(sizePolicy)
            lineEdit_run_p_k7.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k7.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k7.setObjectName("lineEdit_run_p_k7")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k7, 7, 2, 1, 1)
            frame_15 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_15.sizePolicy().hasHeightForWidth())
            frame_15.setSizePolicy(sizePolicy)
            frame_15.setMinimumSize(QtCore.QSize(10, 0))
            frame_15.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_15.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_15.setObjectName("frame_15")
            gridLayout_run_parameters_entry.addWidget(frame_15, 0, 1, 12, 1)
            lineEdit_run_vdna_k10 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k10.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k10.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k10.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k10.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k10.setReadOnly(True)
            lineEdit_run_vdna_k10.setObjectName("lineEdit_run_vdna_k10")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k10, 10, 4, 1, 1)
            label_run_k8 = QtWidgets.QLabel(frame_2)
            label_run_k8.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k8.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k8.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k8.setObjectName("label_run_k8")
            gridLayout_run_parameters_entry.addWidget(label_run_k8, 8, 0, 1, 1)
            lineEdit_run_p_k6 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k6.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k6.setSizePolicy(sizePolicy)
            lineEdit_run_p_k6.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k6.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k6.setObjectName("lineEdit_run_p_k6")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k6, 6, 2, 1, 1)
            label_run_k9 = QtWidgets.QLabel(frame_2)
            label_run_k9.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k9.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k9.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k9.setObjectName("label_run_k9")
            gridLayout_run_parameters_entry.addWidget(label_run_k9, 9, 0, 1, 1)
            frame_10 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_10.sizePolicy().hasHeightForWidth())
            frame_10.setSizePolicy(sizePolicy)
            frame_10.setMinimumSize(QtCore.QSize(10, 0))
            frame_10.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_10.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_10.setObjectName("frame_10")
            gridLayout_run_parameters_entry.addWidget(frame_10, 0, 3, 12, 1)
            lineEdit_run_vdna_k9 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k9.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k9.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k9.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k9.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k9.setReadOnly(True)
            lineEdit_run_vdna_k9.setObjectName("lineEdit_run_vdna_k9")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k9, 9, 4, 1, 1)
            lineEdit_run_vdna_k3 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k3.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k3.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k3.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k3.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k3.setReadOnly(True)
            lineEdit_run_vdna_k3.setObjectName("lineEdit_run_vdna_k3")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k3, 3, 4, 1, 1)
            lineEdit_run_p_k5 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k5.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k5.setSizePolicy(sizePolicy)
            lineEdit_run_p_k5.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k5.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k5.setObjectName("lineEdit_run_p_k5")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k5, 5, 2, 1, 1)
            lineEdit_run_vdna_k2 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k2.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k2.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k2.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k2.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k2.setReadOnly(True)
            lineEdit_run_vdna_k2.setObjectName("lineEdit_run_vdna_k2")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k2, 2, 4, 1, 1)
            label_run_k11 = QtWidgets.QLabel(frame_2)
            label_run_k11.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k11.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k11.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k11.setObjectName("label_run_k11")
            gridLayout_run_parameters_entry.addWidget(label_run_k11, 11, 0, 1, 1)
            lineEdit_run_p_k11 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k11.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k11.setSizePolicy(sizePolicy)
            lineEdit_run_p_k11.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k11.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k11.setObjectName("lineEdit_run_p_k11")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k11, 11, 2, 1, 1)
            label_run_present_values = QtWidgets.QLabel(frame_2)
            label_run_present_values.setAlignment(QtCore.Qt.AlignCenter)
            label_run_present_values.setObjectName("label_run_present_values")
            gridLayout_run_parameters_entry.addWidget(label_run_present_values, 0, 2, 1, 1)
            label_run_k2 = QtWidgets.QLabel(frame_2)
            label_run_k2.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k2.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k2.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k2.setObjectName("label_run_k2")
            gridLayout_run_parameters_entry.addWidget(label_run_k2, 2, 0, 1, 1)
            frame_12 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_12.sizePolicy().hasHeightForWidth())
            frame_12.setSizePolicy(sizePolicy)
            frame_12.setMinimumSize(QtCore.QSize(10, 0))
            frame_12.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_12.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_12.setObjectName("frame_12")
            gridLayout_run_parameters_entry.addWidget(frame_12, 0, 5, 12, 1)
            frame_14 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_14.sizePolicy().hasHeightForWidth())
            frame_14.setSizePolicy(sizePolicy)
            frame_14.setMinimumSize(QtCore.QSize(10, 0))
            frame_14.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_14.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_14.setObjectName("frame_14")
            gridLayout_run_parameters_entry.addWidget(frame_14, 0, 6, 1, 2)
            frame_28 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_28.sizePolicy().hasHeightForWidth())
            frame_28.setSizePolicy(sizePolicy)
            frame_28.setMinimumSize(QtCore.QSize(10, 0))
            frame_28.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_28.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_28.setObjectName("frame_28")
            gridLayout_run_parameters_entry.addWidget(frame_28, 13, 1, 1, 7)
            label_run_k10 = QtWidgets.QLabel(frame_2)
            label_run_k10.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k10.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k10.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k10.setObjectName("label_run_k10")
            gridLayout_run_parameters_entry.addWidget(label_run_k10, 10, 0, 1, 1)
            frame_4 = QtWidgets.QFrame(frame_2)
            frame_4.setMinimumSize(QtCore.QSize(0, 20))
            frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_4.setObjectName("frame_4")
            gridLayout_run_parameters_entry.addWidget(frame_4, 12, 0, 1, 8)
            label_run_k6 = QtWidgets.QLabel(frame_2)
            label_run_k6.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k6.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k6.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k6.setObjectName("label_run_k6")
            gridLayout_run_parameters_entry.addWidget(label_run_k6, 6, 0, 1, 1)
            label_run_k7 = QtWidgets.QLabel(frame_2)
            label_run_k7.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k7.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k7.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k7.setObjectName("label_run_k7")
            gridLayout_run_parameters_entry.addWidget(label_run_k7, 7, 0, 1, 1)
            lineEdit_run_vdna_k4 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k4.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k4.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k4.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k4.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k4.setReadOnly(True)
            lineEdit_run_vdna_k4.setObjectName("lineEdit_run_vdna_k4")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k4, 4, 4, 1, 1)
            label_run_k5 = QtWidgets.QLabel(frame_2)
            label_run_k5.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k5.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k5.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k5.setObjectName("label_run_k5")
            gridLayout_run_parameters_entry.addWidget(label_run_k5, 5, 0, 1, 1)
            lineEdit_run_vdna_k7 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k7.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k7.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k7.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k7.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k7.setReadOnly(True)
            lineEdit_run_vdna_k7.setObjectName("lineEdit_run_vdna_k7")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k7, 7, 4, 1, 1)
            frame_8 = QtWidgets.QFrame(frame_2)
            frame_8.setMinimumSize(QtCore.QSize(0, 20))
            frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_8.setObjectName("frame_8")
            gridLayout_run_parameters_entry.addWidget(frame_8, 15, 0, 1, 8)
            lineEdit_run_temp_p = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_temp_p.sizePolicy().hasHeightForWidth())
            lineEdit_run_temp_p.setSizePolicy(sizePolicy)
            lineEdit_run_temp_p.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_temp_p.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_temp_p.setObjectName("lineEdit_run_temp_p")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_temp_p, 14, 2, 1, 1)
            frame_26 = QtWidgets.QFrame(frame_2)
            frame_26.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_26.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_26.setObjectName("frame_26")
            gridLayout_run_parameters_entry.addWidget(frame_26, 14, 5, 1, 1)
            frame_5 = QtWidgets.QFrame(frame_2)
            frame_5.setMinimumSize(QtCore.QSize(0, 20))
            frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_5.setObjectName("frame_5")
            gridLayout_run_parameters_entry.addWidget(frame_5, 19, 0, 1, 8)
            label_run_temp = QtWidgets.QLabel(frame_2)
            label_run_temp.setMinimumSize(QtCore.QSize(0, 25))
            label_run_temp.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_temp.setAlignment(QtCore.Qt.AlignCenter)
            label_run_temp.setObjectName("label_run_temp")
            gridLayout_run_parameters_entry.addWidget(label_run_temp, 14, 0, 1, 1)
            frame_23 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_23.sizePolicy().hasHeightForWidth())
            frame_23.setSizePolicy(sizePolicy)
            frame_23.setMinimumSize(QtCore.QSize(10, 0))
            frame_23.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_23.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_23.setObjectName("frame_23")
            gridLayout_run_parameters_entry.addWidget(frame_23, 14, 1, 1, 1)
            frame_27 = QtWidgets.QFrame(frame_2)
            frame_27.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_27.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_27.setObjectName("frame_27")
            gridLayout_run_parameters_entry.addWidget(frame_27, 14, 6, 1, 2)
            lineEdit_run_temp_vdna = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_temp_vdna.sizePolicy().hasHeightForWidth())
            lineEdit_run_temp_vdna.setSizePolicy(sizePolicy)
            lineEdit_run_temp_vdna.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_temp_vdna.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_temp_vdna.setReadOnly(True)
            lineEdit_run_temp_vdna.setObjectName("lineEdit_run_temp_vdna")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_temp_vdna, 14, 4, 1, 1)
            label_run_tempereture = QtWidgets.QLabel(frame_2)
            label_run_tempereture.setAlignment(QtCore.Qt.AlignCenter)
            label_run_tempereture.setObjectName("label_run_tempereture")
            gridLayout_run_parameters_entry.addWidget(label_run_tempereture, 13, 0, 1, 1)
            label_run_kinetic_values = QtWidgets.QLabel(frame_2)
            label_run_kinetic_values.setAlignment(QtCore.Qt.AlignCenter)
            label_run_kinetic_values.setObjectName("label_run_kinetic_values")
            gridLayout_run_parameters_entry.addWidget(label_run_kinetic_values, 0, 0, 1, 1)
            pushButton_7 = QtWidgets.QPushButton(frame_2)
            pushButton_7.setObjectName("pushButton_7")
            gridLayout_run_parameters_entry.addWidget(pushButton_7, 8, 6, 1, 2)
            pushButton_reset_to_originals = QtWidgets.QPushButton(frame_2)
            pushButton_reset_to_originals.setObjectName("pushButton_reset_to_originals")
            gridLayout_run_parameters_entry.addWidget(pushButton_reset_to_originals, 7, 6, 1, 2)
            pushButton_set_vdna_defaults = QtWidgets.QPushButton(frame_2)
            pushButton_set_vdna_defaults.setObjectName("pushButton_set_vdna_defaults")
            gridLayout_run_parameters_entry.addWidget(pushButton_set_vdna_defaults, 6, 6, 1, 2)
            frame_25 = QtWidgets.QFrame(frame_2)
            frame_25.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_25.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_25.setObjectName("frame_25")
            gridLayout_run_parameters_entry.addWidget(frame_25, 14, 3, 1, 1)
            lineEdit_run_p_k3 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k3.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k3.setSizePolicy(sizePolicy)
            lineEdit_run_p_k3.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k3.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k3.setObjectName("lineEdit_run_p_k3")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k3, 3, 2, 1, 1)
            label_run_k1 = QtWidgets.QLabel(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(label_run_k1.sizePolicy().hasHeightForWidth())
            label_run_k1.setSizePolicy(sizePolicy)
            label_run_k1.setMinimumSize(QtCore.QSize(0, 25))
            label_run_k1.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_run_k1.setAlignment(QtCore.Qt.AlignCenter)
            label_run_k1.setObjectName("label_run_k1")
            gridLayout_run_parameters_entry.addWidget(label_run_k1, 1, 0, 1, 1)
            lineEdit_run_p_k8 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_p_k8.sizePolicy().hasHeightForWidth())
            lineEdit_run_p_k8.setSizePolicy(sizePolicy)
            lineEdit_run_p_k8.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_p_k8.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_p_k8.setObjectName("lineEdit_run_p_k8")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_p_k8, 8, 2, 1, 1)
            lineEdit_run_vdna_k11 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k11.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k11.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k11.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k11.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k11.setReadOnly(True)
            lineEdit_run_vdna_k11.setObjectName("lineEdit_run_vdna_k11")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k11, 11, 4, 1, 1)
            lineEdit_run_vdna_k6 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_vdna_k6.sizePolicy().hasHeightForWidth())
            lineEdit_run_vdna_k6.setSizePolicy(sizePolicy)
            lineEdit_run_vdna_k6.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_vdna_k6.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_vdna_k6.setReadOnly(True)
            lineEdit_run_vdna_k6.setObjectName("lineEdit_run_vdna_k6")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_vdna_k6, 6, 4, 1, 1)
            frame_30 = QtWidgets.QFrame(frame_2)
            frame_30.setMinimumSize(QtCore.QSize(0, 20))
            frame_30.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_30.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_30.setObjectName("frame_30")
            gridLayout_run_parameters_entry.addWidget(frame_30, 16, 2, 1, 1)
            frame_32 = QtWidgets.QFrame(frame_2)
            frame_32.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_32.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_32.setObjectName("frame_32")
            gridLayout_run_parameters_entry.addWidget(frame_32, 16, 6, 1, 1)
            frame_31 = QtWidgets.QFrame(frame_2)
            frame_31.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_31.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_31.setObjectName("frame_31")
            gridLayout_run_parameters_entry.addWidget(frame_31, 16, 4, 1, 1)
            frame_21 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_21.sizePolicy().hasHeightForWidth())
            frame_21.setSizePolicy(sizePolicy)
            frame_21.setMinimumSize(QtCore.QSize(10, 0))
            frame_21.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_21.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_21.setObjectName("frame_21")
            gridLayout_run_parameters_entry.addWidget(frame_21, 16, 1, 3, 1)
            frame_19 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_19.sizePolicy().hasHeightForWidth())
            frame_19.setSizePolicy(sizePolicy)
            frame_19.setMinimumSize(QtCore.QSize(10, 0))
            frame_19.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_19.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_19.setObjectName("frame_19")
            gridLayout_run_parameters_entry.addWidget(frame_19, 16, 3, 3, 1)
            frame_17 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_17.sizePolicy().hasHeightForWidth())
            frame_17.setSizePolicy(sizePolicy)
            frame_17.setMinimumSize(QtCore.QSize(10, 0))
            frame_17.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_17.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_17.setObjectName("frame_17")
            gridLayout_run_parameters_entry.addWidget(frame_17, 16, 5, 3, 1)
            label_run_d_temp = QtWidgets.QLabel(frame_2)
            label_run_d_temp.setAlignment(QtCore.Qt.AlignCenter)
            label_run_d_temp.setObjectName("label_run_d_temp")
            gridLayout_run_parameters_entry.addWidget(label_run_d_temp, 17, 2, 1, 1)
            lineEdit_run_start_temp = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_start_temp.sizePolicy().hasHeightForWidth())
            lineEdit_run_start_temp.setSizePolicy(sizePolicy)
            lineEdit_run_start_temp.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_start_temp.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_start_temp.setObjectName("lineEdit_run_start_temp")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_start_temp, 18, 0, 1, 1)
            lineEdit_run_d_temp = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_d_temp.sizePolicy().hasHeightForWidth())
            lineEdit_run_d_temp.setSizePolicy(sizePolicy)
            lineEdit_run_d_temp.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_d_temp.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_d_temp.setObjectName("lineEdit_run_d_temp")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_d_temp, 18, 2, 1, 1)
            label_run_end_temp = QtWidgets.QLabel(frame_2)
            label_run_end_temp.setAlignment(QtCore.Qt.AlignCenter)
            label_run_end_temp.setObjectName("label_run_end_temp")
            gridLayout_run_parameters_entry.addWidget(label_run_end_temp, 17, 4, 1, 1)
            lineEdit_run_end_temp = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_end_temp.sizePolicy().hasHeightForWidth())
            lineEdit_run_end_temp.setSizePolicy(sizePolicy)
            lineEdit_run_end_temp.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_end_temp.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_end_temp.setObjectName("lineEdit_run_end_temp")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_end_temp, 18, 4, 1, 1)
            label_run_n_stages = QtWidgets.QLabel(frame_2)
            label_run_n_stages.setAlignment(QtCore.Qt.AlignCenter)
            label_run_n_stages.setObjectName("label_run_n_stages")
            gridLayout_run_parameters_entry.addWidget(label_run_n_stages, 17, 6, 1, 1)
            lineEdit_run_total_time = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_total_time.sizePolicy().hasHeightForWidth())
            lineEdit_run_total_time.setSizePolicy(sizePolicy)
            lineEdit_run_total_time.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_total_time.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_total_time.setReadOnly(True)
            lineEdit_run_total_time.setObjectName("lineEdit_run_total_time")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_total_time, 21, 6, 1, 1)
            label_run_dumps_per_stage = QtWidgets.QLabel(frame_2)
            label_run_dumps_per_stage.setAlignment(QtCore.Qt.AlignCenter)
            label_run_dumps_per_stage.setObjectName("label_run_dumps_per_stage")
            gridLayout_run_parameters_entry.addWidget(label_run_dumps_per_stage, 20, 4, 1, 1)
            label_run_start_time = QtWidgets.QLabel(frame_2)
            label_run_start_time.setAlignment(QtCore.Qt.AlignCenter)
            label_run_start_time.setObjectName("label_run_start_time")
            gridLayout_run_parameters_entry.addWidget(label_run_start_time, 20, 0, 1, 1)
            label_run_total_time = QtWidgets.QLabel(frame_2)
            label_run_total_time.setAlignment(QtCore.Qt.AlignCenter)
            label_run_total_time.setObjectName("label_run_total_time")
            gridLayout_run_parameters_entry.addWidget(label_run_total_time, 20, 6, 1, 1)
            lineEdit_run_dumps_per_stage = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_dumps_per_stage.sizePolicy().hasHeightForWidth())
            lineEdit_run_dumps_per_stage.setSizePolicy(sizePolicy)
            lineEdit_run_dumps_per_stage.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_dumps_per_stage.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_dumps_per_stage.setObjectName("lineEdit_run_dumps_per_stage")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_dumps_per_stage, 21, 4, 1, 1)
            label_run_time_per_stage = QtWidgets.QLabel(frame_2)
            label_run_time_per_stage.setAlignment(QtCore.Qt.AlignCenter)
            label_run_time_per_stage.setObjectName("label_run_time_per_stage")
            gridLayout_run_parameters_entry.addWidget(label_run_time_per_stage, 20, 2, 1, 1)
            lineEdit_run_time_per_stage = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_run_time_per_stage.sizePolicy().hasHeightForWidth())
            lineEdit_run_time_per_stage.setSizePolicy(sizePolicy)
            lineEdit_run_time_per_stage.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_run_time_per_stage.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_run_time_per_stage.setObjectName("lineEdit_run_time_per_stage")
            gridLayout_run_parameters_entry.addWidget(lineEdit_run_time_per_stage, 21, 2, 1, 1)
            label_run_start_temp = QtWidgets.QLabel(frame_2)
            label_run_start_temp.setAlignment(QtCore.Qt.AlignCenter)
            label_run_start_temp.setObjectName("label_run_start_temp")
            gridLayout_run_parameters_entry.addWidget(label_run_start_temp, 17, 0, 1, 1)
            checkBox_annealing = QtWidgets.QCheckBox(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(checkBox_annealing.sizePolicy().hasHeightForWidth())
            checkBox_annealing.setSizePolicy(sizePolicy)
            checkBox_annealing.setMinimumSize(QtCore.QSize(100, 0))
            font = QtGui.QFont()
            font.setBold(True)
            font.setWeight(75)
            checkBox_annealing.setFont(font)
            checkBox_annealing.setObjectName("checkBox_annealing")
            gridLayout_run_parameters_entry.addWidget(checkBox_annealing, 16, 0, 1, 1)
            frame_35 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_35.sizePolicy().hasHeightForWidth())
            frame_35.setSizePolicy(sizePolicy)
            frame_35.setMinimumSize(QtCore.QSize(10, 0))
            frame_35.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_35.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_35.setObjectName("frame_35")
            gridLayout_run_parameters_entry.addWidget(frame_35, 20, 3, 3, 1)
            frame_37 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_37.sizePolicy().hasHeightForWidth())
            frame_37.setSizePolicy(sizePolicy)
            frame_37.setMinimumSize(QtCore.QSize(10, 0))
            frame_37.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_37.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_37.setObjectName("frame_37")
            gridLayout_run_parameters_entry.addWidget(frame_37, 20, 5, 3, 1)
            frame_33 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_33.sizePolicy().hasHeightForWidth())
            frame_33.setSizePolicy(sizePolicy)
            frame_33.setMinimumSize(QtCore.QSize(10, 0))
            frame_33.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_33.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_33.setObjectName("frame_33")
            gridLayout_run_parameters_entry.addWidget(frame_33, 20, 1, 3, 1)
            verticalLayout.addLayout(gridLayout_run_parameters_entry)
            frame = QtWidgets.QFrame(frame_2)
            frame.setMinimumSize(QtCore.QSize(0, 20))
            frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame.setFrameShadow(QtWidgets.QFrame.Raised)
            frame.setObjectName("frame")
            verticalLayout.addWidget(frame)
            horizontalLayout_save_cancel = QtWidgets.QHBoxLayout()
            horizontalLayout_save_cancel.setObjectName("horizontalLayout_save_cancel")
            buttonBox_write_edit = QtWidgets.QDialogButtonBox(frame_2)
            font = QtGui.QFont()
            font.setKerning(True)
            buttonBox_write_edit.setFont(font)
            buttonBox_write_edit.setOrientation(QtCore.Qt.Horizontal)
            buttonBox_write_edit.setStandardButtons(QtWidgets.QDialogButtonBox.Save)
            buttonBox_write_edit.setObjectName("buttonBox_write_edit")
            horizontalLayout_save_cancel.addWidget(buttonBox_write_edit)
            pushButton_cancel = QtWidgets.QPushButton(frame_2)
            pushButton_cancel.setObjectName("pushButton_cancel")
            horizontalLayout_save_cancel.addWidget(pushButton_cancel)
            verticalLayout.addLayout(horizontalLayout_save_cancel)
            gridLayout.addLayout(verticalLayout, 0, 0, 1, 1)
            gridLayout_2.addWidget(frame_2, 0, 0, 1, 1)

            def retranslateUi(Dialog):
                _translate = QtCore.QCoreApplication.translate
                Dialog.setWindowTitle(_translate("Dialog", "Advanced run parameters"))
                label_run_vdna_defaults.setText(_translate("Dialog", "VDNA defaults"))
                label_run_k4.setText(_translate("Dialog", "k 4"))
                label_run_k3.setText(_translate("Dialog", "k 3"))
                label_run_k8.setText(_translate("Dialog", "k 8"))
                label_run_k9.setText(_translate("Dialog", "k 9"))
                label_run_k11.setText(_translate("Dialog", "k 11"))
                label_run_present_values.setText(_translate("Dialog", "Present values"))
                label_run_k2.setText(_translate("Dialog", "k 2"))
                label_run_k10.setText(_translate("Dialog", "k 10"))
                label_run_k6.setText(_translate("Dialog", "k 6"))
                label_run_k7.setText(_translate("Dialog", "k 7"))
                label_run_k5.setText(_translate("Dialog", "k 5"))
                label_run_temp.setText(_translate("Dialog", "Temp"))
                label_run_tempereture.setText(_translate("Dialog", "Temperature"))
                label_run_kinetic_values.setText(_translate("Dialog", "Kinetic values"))
                pushButton_7.setText(_translate("Dialog", "????"))
                pushButton_reset_to_originals.setText(_translate("Dialog", "Reset to originals"))
                pushButton_set_vdna_defaults.setText(_translate("Dialog", "Set VDNA defaults"))
                label_run_k1.setText(_translate("Dialog", "k 1"))
                label_run_d_temp.setText(_translate("Dialog", "Δ temperature"))
                label_run_end_temp.setText(_translate("Dialog", "End temperature"))
                label_run_n_stages.setText(_translate("Dialog", "# stages"))
                label_run_dumps_per_stage.setText(_translate("Dialog", "# dumps per stage"))
                label_run_start_time.setText(_translate("Dialog", "Start time"))
                label_run_total_time.setText(_translate("Dialog", "Total time"))
                label_run_time_per_stage.setText(_translate("Dialog", "Time per stage"))
                label_run_start_temp.setText(_translate("Dialog", "Start temperature"))
                checkBox_annealing.setText(_translate("Dialog", "Annealing"))
                pushButton_cancel.setText(_translate("Dialog", "Cancel"))

            def close_dialog():
                run_advanced_dialog.reject()

            retranslateUi(run_advanced_dialog)
            buttonBox_write_edit.rejected.connect(run_advanced_dialog.reject)
            buttonBox_write_edit.accepted.connect(run_advanced_dialog.accept)

            kt_n_ann_cells = {'Temp_p': lineEdit_run_temp_p,
                              'k1_p': lineEdit_run_p_k1,
                              'k2_P': lineEdit_run_p_k2,
                              'k3_p': lineEdit_run_p_k3,
                              'k4_p': lineEdit_run_p_k4,
                              'k5_p': lineEdit_run_p_k5,
                              'k6_p': lineEdit_run_p_k6,
                              'k7_p': lineEdit_run_p_k7,
                              'k8_p': lineEdit_run_p_k8,
                              'k9_P': lineEdit_run_p_k9,
                              'k10_p': lineEdit_run_p_k10,
                              'k11_p': lineEdit_run_p_k11,
                              'run_start_time': lineEdit_run_start_time,
                              'run_time_per_stage': lineEdit_run_time_per_stage,
                              'run_dumps_per_stage': lineEdit_run_dumps_per_stage,
                              'run_start_temp': lineEdit_run_start_temp,
                              'run_d_temp': lineEdit_run_d_temp,
                              'run_end_temp': lineEdit_run_end_temp}

            set_k_values(source_file, 'get')
            pushButton_set_vdna_defaults.clicked.connect(lambda: set_k_values(source_file, 'set_d'))
            pushButton_reset_to_originals.clicked.connect(lambda: set_k_values(source_file, 'set_o'))
            pushButton_cancel.clicked.connect(close_dialog)
            checkBox_annealing.clicked.connect(lambda: activate_annealing('u_clicked'))

            lineEdit_run_temp_p.textChanged.connect(lambda: validate_cell_values('Temp_p'))
            lineEdit_run_p_k1.textChanged.connect(lambda: validate_cell_values('k1_p'))
            lineEdit_run_p_k2.textChanged.connect(lambda: validate_cell_values('k2_P'))
            lineEdit_run_p_k3.textChanged.connect(lambda: validate_cell_values('k3_p'))
            lineEdit_run_p_k4.textChanged.connect(lambda: validate_cell_values('k4_p'))
            lineEdit_run_p_k5.textChanged.connect(lambda: validate_cell_values('k5_p'))
            lineEdit_run_p_k6.textChanged.connect(lambda: validate_cell_values('k6_p'))
            lineEdit_run_p_k7.textChanged.connect(lambda: validate_cell_values('k7_p'))
            lineEdit_run_p_k8.textChanged.connect(lambda: validate_cell_values('k8_p'))
            lineEdit_run_p_k9.textChanged.connect(lambda: validate_cell_values('k9_P'))
            lineEdit_run_p_k10.textChanged.connect(lambda: validate_cell_values('k10_p'))
            lineEdit_run_p_k11.textChanged.connect(lambda: validate_cell_values('k11_p'))
            lineEdit_run_start_time.textChanged.connect(lambda: validate_cell_values('run_start_time'))
            lineEdit_run_time_per_stage.textChanged.connect(lambda: validate_cell_values('run_time_per_stage'))
            lineEdit_run_dumps_per_stage.textChanged.connect(lambda: validate_cell_values('run_dumps_per_stage'))
            lineEdit_run_start_temp.textChanged.connect(lambda: validate_cell_values('run_start_temp'))
            lineEdit_run_d_temp.textChanged.connect(lambda: validate_cell_values('run_d_temp'))
            lineEdit_run_end_temp.textChanged.connect(lambda: validate_cell_values('run_end_temp'))

            activate_annealing('')
            set_previous()

            run_run_advanced_dialog = run_advanced_dialog.exec_()

            if run_run_advanced_dialog == 1:
                global run_annealing_g, run_kt_g

                if checkBox_annealing.isChecked():
                    run_annealing_g = True
                elif not checkBox_annealing.isChecked():
                    run_annealing_g = False

                run_kt_g = True
                save_rnf_data()
                self.checkBox_run_advanced.setChecked(True)
                run_activate_run_advanced()

        def write_activate_write_advanced():
            if self.checkBox_write_advanced.isChecked():
                if write_kt_status_g:
                    self.label_messages.setText('')

                elif not write_kt_status_g:
                    self.label_messages.setText('Advanced parameters not entered in advanced menu ✘')
                    self.checkBox_write_advanced.setChecked(False)

        def write_advanced():

            source_file = self.lineEdit_write_browse_source.text()

            def update_save_button():

                all_cells_check = sum([1 for i in kt_cells.values() if i.text().endswith('.')]) == 0

                if all_cells_check:
                    buttonBox_write_edit.setEnabled(True)

                elif not all_cells_check:
                    buttonBox_write_edit.setDisabled(True)

            def set_k_values(f_path, call_from):

                cells = {'Temp': [lineEdit_write_temp_p, lineEdit_write_temp_vdna],
                         'k1': [lineEdit_write_p_k1, lineEdit_write_vdna_k1],
                         'k2': [lineEdit_write_p_k2, lineEdit_write_vdna_k2],
                         'k3': [lineEdit_write_p_k3, lineEdit_write_vdna_k3],
                         'k4': [lineEdit_write_p_k4, lineEdit_write_vdna_k4],
                         'k5': [lineEdit_write_p_k5, lineEdit_write_vdna_k5],
                         'k6': [lineEdit_write_p_k6, lineEdit_write_vdna_k6],
                         'k7': [lineEdit_write_p_k7, lineEdit_write_vdna_k7],
                         'k8': [lineEdit_write_p_k8, lineEdit_write_vdna_k8],
                         'k9': [lineEdit_write_p_k9, lineEdit_write_vdna_k9],
                         'k10': [lineEdit_write_p_k10, lineEdit_write_vdna_k10],
                         'k11': [lineEdit_write_p_k11, lineEdit_write_vdna_k11]}

                def set_file_write_parameters(a):
                    k_data_r = get_k_values('system_files/reference_files/ref_file.bngl')

                    if a == 'get':
                        k_data_p = get_k_values(f_path) if write_kt_list_g == {} else write_kt_list_g
                        for ii in cells.items():
                            if ii[0] in k_data_p:
                                ii[1][0].setText(k_data_p[ii[0]])
                            else:
                                ii[1][0].setText('')

                            if ii[0] in k_data_r:
                                ii[1][1].setText(k_data_r[ii[0]])
                            else:
                                ii[1][1].setText('')

                    elif a == 'set_d':
                        for ii in cells.items():
                            if ii[0] in k_data_r:
                                ii[1][0].setText(k_data_r[ii[0]])
                            else:
                                ii[1][0].setText('')

                    elif a == 'set_o':
                        k_data_p = get_k_values(f_path)
                        for ii in cells.items():
                            if ii[0] in k_data_p:
                                ii[1][0].setText(k_data_p[ii[0]])
                            else:
                                ii[1][0].setText('')

                set_file_write_parameters(call_from)

            def save_kt_data():
                global write_kt_list_g

                k_list_values = {'Temp': lineEdit_write_temp_p.text(),
                                 'k1': lineEdit_write_p_k1.text(), 'k2': lineEdit_write_p_k2.text(),
                                 'k3': lineEdit_write_p_k3.text(), 'k4': lineEdit_write_p_k4.text(),
                                 'k5': lineEdit_write_p_k5.text(), 'k6': lineEdit_write_p_k6.text(),
                                 'k7': lineEdit_write_p_k7.text(), 'k8': lineEdit_write_p_k8.text(),
                                 'k9': lineEdit_write_p_k9.text(), 'k10': lineEdit_write_p_k10.text(),
                                 'k11': lineEdit_write_p_k11.text()}

                del_k = []
                for key, val in k_list_values.items():
                    if val == '':
                        del_k.append(key)

                for k in del_k:
                    del k_list_values[k]

                write_kt_list_g = k_list_values

            def validate_cell_values(call_from):
                global write_n_run_kt_validated_g

                current_item = kt_cells[call_from].text()

                write_n_run_kt_validated_g = current_item
                validate_num(write_n_run_kt_validate, '')

                if current_item != write_n_run_kt_validated_g:
                    kt_cells[call_from].setText(write_n_run_kt_validated_g)
                    write_n_run_kt_validated_g = ''

                update_save_button()

            write_advanced_dialog = QDialog()
            write_advanced_dialog.setObjectName("Dialog")
            write_advanced_dialog.resize(600, 610)
            write_advanced_dialog.setMaximumSize(QtCore.QSize(700, 810))
            gridLayout_2 = QtWidgets.QGridLayout(write_advanced_dialog)
            gridLayout_2.setObjectName("gridLayout_2")
            frame_2 = QtWidgets.QFrame(write_advanced_dialog)
            frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_2.setObjectName("frame_2")
            gridLayout = QtWidgets.QGridLayout(frame_2)
            gridLayout.setContentsMargins(20, 20, 20, 20)
            gridLayout.setObjectName("gridLayout")
            verticalLayout = QtWidgets.QVBoxLayout()
            verticalLayout.setObjectName("verticalLayout")
            gridLayout_write_parameters_entry = QtWidgets.QGridLayout()
            gridLayout_write_parameters_entry.setObjectName("gridLayout_write_parameters_entry")
            label_write_vdna_defaults = QtWidgets.QLabel(frame_2)
            label_write_vdna_defaults.setAlignment(QtCore.Qt.AlignCenter)
            label_write_vdna_defaults.setObjectName("label_write_vdna_defaults")
            gridLayout_write_parameters_entry.addWidget(label_write_vdna_defaults, 0, 4, 1, 1)
            lineEdit_write_p_k4 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k4.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k4.setSizePolicy(sizePolicy)
            lineEdit_write_p_k4.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k4.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k4.setObjectName("lineEdit_write_p_k4")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k4, 4, 2, 1, 1)
            lineEdit_write_vdna_k1 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k1.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k1.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k1.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k1.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k1.setReadOnly(True)
            lineEdit_write_vdna_k1.setObjectName("lineEdit_write_vdna_k1")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k1, 1, 4, 1, 1)
            lineEdit_write_p_k10 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k10.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k10.setSizePolicy(sizePolicy)
            lineEdit_write_p_k10.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k10.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k10.setObjectName("lineEdit_write_p_k10")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k10, 10, 2, 1, 1)
            lineEdit_write_p_k9 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k9.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k9.setSizePolicy(sizePolicy)
            lineEdit_write_p_k9.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k9.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k9.setObjectName("lineEdit_write_p_k9")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k9, 9, 2, 1, 1)
            label_write_k4 = QtWidgets.QLabel(frame_2)
            label_write_k4.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k4.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k4.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k4.setObjectName("label_write_k4")
            gridLayout_write_parameters_entry.addWidget(label_write_k4, 4, 0, 1, 1)
            lineEdit_write_vdna_k5 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k5.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k5.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k5.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k5.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k5.setReadOnly(True)
            lineEdit_write_vdna_k5.setObjectName("lineEdit_write_vdna_k5")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k5, 5, 4, 1, 1)
            lineEdit_write_vdna_k8 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k8.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k8.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k8.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k8.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k8.setReadOnly(True)
            lineEdit_write_vdna_k8.setObjectName("lineEdit_write_vdna_k8")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k8, 8, 4, 1, 1)
            label_write_k8 = QtWidgets.QLabel(frame_2)
            label_write_k8.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k8.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k8.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k8.setObjectName("label_write_k8")
            gridLayout_write_parameters_entry.addWidget(label_write_k8, 8, 0, 1, 1)
            lineEdit_write_vdna_k9 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k9.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k9.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k9.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k9.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k9.setReadOnly(True)
            lineEdit_write_vdna_k9.setObjectName("lineEdit_write_vdna_k9")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k9, 9, 4, 1, 1)
            lineEdit_write_vdna_k2 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k2.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k2.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k2.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k2.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k2.setReadOnly(True)
            lineEdit_write_vdna_k2.setObjectName("lineEdit_write_vdna_k2")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k2, 2, 4, 1, 1)
            lineEdit_write_p_k6 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k6.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k6.setSizePolicy(sizePolicy)
            lineEdit_write_p_k6.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k6.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k6.setObjectName("lineEdit_write_p_k6")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k6, 6, 2, 1, 1)
            lineEdit_write_p_k2 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k2.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k2.setSizePolicy(sizePolicy)
            lineEdit_write_p_k2.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k2.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k2.setObjectName("lineEdit_write_p_k2")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k2, 2, 2, 1, 1)
            lineEdit_write_vdna_k3 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k3.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k3.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k3.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k3.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k3.setReadOnly(True)
            lineEdit_write_vdna_k3.setObjectName("lineEdit_write_vdna_k3")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k3, 3, 4, 1, 1)
            lineEdit_write_p_k1 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k1.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k1.setSizePolicy(sizePolicy)
            lineEdit_write_p_k1.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k1.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k1.setObjectName("lineEdit_write_p_k1")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k1, 1, 2, 1, 1)
            label_write_k3 = QtWidgets.QLabel(frame_2)
            label_write_k3.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k3.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k3.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k3.setObjectName("label_write_k3")
            gridLayout_write_parameters_entry.addWidget(label_write_k3, 3, 0, 1, 1)
            lineEdit_write_p_k5 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k5.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k5.setSizePolicy(sizePolicy)
            lineEdit_write_p_k5.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k5.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k5.setObjectName("lineEdit_write_p_k5")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k5, 5, 2, 1, 1)
            lineEdit_write_p_k7 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k7.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k7.setSizePolicy(sizePolicy)
            lineEdit_write_p_k7.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k7.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k7.setObjectName("lineEdit_write_p_k7")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k7, 7, 2, 1, 1)
            frame_15 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_15.sizePolicy().hasHeightForWidth())
            frame_15.setSizePolicy(sizePolicy)
            frame_15.setMinimumSize(QtCore.QSize(10, 0))
            frame_15.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_15.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_15.setObjectName("frame_15")
            gridLayout_write_parameters_entry.addWidget(frame_15, 0, 1, 12, 1)
            lineEdit_write_vdna_k10 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k10.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k10.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k10.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k10.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k10.setReadOnly(True)
            lineEdit_write_vdna_k10.setObjectName("lineEdit_write_vdna_k10")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k10, 10, 4, 1, 1)
            label_write_k9 = QtWidgets.QLabel(frame_2)
            label_write_k9.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k9.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k9.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k9.setObjectName("label_write_k9")
            gridLayout_write_parameters_entry.addWidget(label_write_k9, 9, 0, 1, 1)
            label_write_k11 = QtWidgets.QLabel(frame_2)
            label_write_k11.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k11.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k11.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k11.setObjectName("label_write_k11")
            gridLayout_write_parameters_entry.addWidget(label_write_k11, 11, 0, 1, 1)
            frame_10 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_10.sizePolicy().hasHeightForWidth())
            frame_10.setSizePolicy(sizePolicy)
            frame_10.setMinimumSize(QtCore.QSize(10, 0))
            frame_10.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_10.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_10.setObjectName("frame_10")
            gridLayout_write_parameters_entry.addWidget(frame_10, 0, 3, 12, 1)
            lineEdit_write_p_k11 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k11.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k11.setSizePolicy(sizePolicy)
            lineEdit_write_p_k11.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k11.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k11.setObjectName("lineEdit_write_p_k11")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k11, 11, 2, 1, 1)
            label_write_k2 = QtWidgets.QLabel(frame_2)
            label_write_k2.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k2.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k2.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k2.setObjectName("label_write_k2")
            gridLayout_write_parameters_entry.addWidget(label_write_k2, 2, 0, 1, 1)
            label_write_present_values = QtWidgets.QLabel(frame_2)
            label_write_present_values.setAlignment(QtCore.Qt.AlignCenter)
            label_write_present_values.setObjectName("label_write_present_values")
            gridLayout_write_parameters_entry.addWidget(label_write_present_values, 0, 2, 1, 1)
            label_write_tempereture = QtWidgets.QLabel(frame_2)
            label_write_tempereture.setAlignment(QtCore.Qt.AlignCenter)
            label_write_tempereture.setObjectName("label_write_tempereture")
            gridLayout_write_parameters_entry.addWidget(label_write_tempereture, 13, 0, 1, 1)
            pushButton_set_vdna_defaults = QtWidgets.QPushButton(frame_2)
            pushButton_set_vdna_defaults.setObjectName("pushButton_set_vdna_defaults")
            gridLayout_write_parameters_entry.addWidget(pushButton_set_vdna_defaults, 6, 6, 1, 2)
            lineEdit_write_vdna_k6 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k6.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k6.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k6.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k6.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k6.setReadOnly(True)
            lineEdit_write_vdna_k6.setObjectName("lineEdit_write_vdna_k6")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k6, 6, 4, 1, 1)
            pushButton_7 = QtWidgets.QPushButton(frame_2)
            pushButton_7.setObjectName("pushButton_7")
            gridLayout_write_parameters_entry.addWidget(pushButton_7, 8, 6, 1, 2)
            label_write_k1 = QtWidgets.QLabel(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(label_write_k1.sizePolicy().hasHeightForWidth())
            label_write_k1.setSizePolicy(sizePolicy)
            label_write_k1.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k1.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k1.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k1.setObjectName("label_write_k1")
            gridLayout_write_parameters_entry.addWidget(label_write_k1, 1, 0, 1, 1)
            lineEdit_write_p_k8 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k8.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k8.setSizePolicy(sizePolicy)
            lineEdit_write_p_k8.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k8.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k8.setObjectName("lineEdit_write_p_k8")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k8, 8, 2, 1, 1)
            lineEdit_write_temp_vdna = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_temp_vdna.sizePolicy().hasHeightForWidth())
            lineEdit_write_temp_vdna.setSizePolicy(sizePolicy)
            lineEdit_write_temp_vdna.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_temp_vdna.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_temp_vdna.setReadOnly(True)
            lineEdit_write_temp_vdna.setObjectName("lineEdit_write_temp_vdna")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_temp_vdna, 14, 4, 1, 1)
            lineEdit_write_vdna_k11 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k11.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k11.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k11.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k11.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k11.setReadOnly(True)
            lineEdit_write_vdna_k11.setObjectName("lineEdit_write_vdna_k11")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k11, 11, 4, 1, 1)
            pushButton_reset_to_originals = QtWidgets.QPushButton(frame_2)
            pushButton_reset_to_originals.setObjectName("pushButton_reset_to_originals")
            gridLayout_write_parameters_entry.addWidget(pushButton_reset_to_originals, 7, 6, 1, 2)
            label_write_kinetic_values = QtWidgets.QLabel(frame_2)
            label_write_kinetic_values.setAlignment(QtCore.Qt.AlignCenter)
            label_write_kinetic_values.setObjectName("label_write_kinetic_values")
            gridLayout_write_parameters_entry.addWidget(label_write_kinetic_values, 0, 0, 1, 1)
            frame_25 = QtWidgets.QFrame(frame_2)
            frame_25.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_25.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_25.setObjectName("frame_25")
            gridLayout_write_parameters_entry.addWidget(frame_25, 14, 3, 1, 1)
            lineEdit_write_p_k3 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_p_k3.sizePolicy().hasHeightForWidth())
            lineEdit_write_p_k3.setSizePolicy(sizePolicy)
            lineEdit_write_p_k3.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_p_k3.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_p_k3.setObjectName("lineEdit_write_p_k3")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_p_k3, 3, 2, 1, 1)
            lineEdit_write_vdna_k4 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k4.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k4.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k4.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k4.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k4.setReadOnly(True)
            lineEdit_write_vdna_k4.setObjectName("lineEdit_write_vdna_k4")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k4, 4, 4, 1, 1)
            frame_12 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_12.sizePolicy().hasHeightForWidth())
            frame_12.setSizePolicy(sizePolicy)
            frame_12.setMinimumSize(QtCore.QSize(10, 0))
            frame_12.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_12.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_12.setObjectName("frame_12")
            gridLayout_write_parameters_entry.addWidget(frame_12, 0, 5, 12, 1)
            frame_4 = QtWidgets.QFrame(frame_2)
            frame_4.setMinimumSize(QtCore.QSize(0, 20))
            frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_4.setObjectName("frame_4")
            gridLayout_write_parameters_entry.addWidget(frame_4, 12, 0, 1, 8)
            lineEdit_write_vdna_k7 = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_vdna_k7.sizePolicy().hasHeightForWidth())
            lineEdit_write_vdna_k7.setSizePolicy(sizePolicy)
            lineEdit_write_vdna_k7.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_vdna_k7.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_vdna_k7.setReadOnly(True)
            lineEdit_write_vdna_k7.setObjectName("lineEdit_write_vdna_k7")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_vdna_k7, 7, 4, 1, 1)
            label_write_k10 = QtWidgets.QLabel(frame_2)
            label_write_k10.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k10.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k10.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k10.setObjectName("label_write_k10")
            gridLayout_write_parameters_entry.addWidget(label_write_k10, 10, 0, 1, 1)
            lineEdit_write_temp_p = QtWidgets.QLineEdit(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(lineEdit_write_temp_p.sizePolicy().hasHeightForWidth())
            lineEdit_write_temp_p.setSizePolicy(sizePolicy)
            lineEdit_write_temp_p.setMinimumSize(QtCore.QSize(100, 0))
            lineEdit_write_temp_p.setAlignment(QtCore.Qt.AlignCenter)
            lineEdit_write_temp_p.setObjectName("lineEdit_write_temp_p")
            gridLayout_write_parameters_entry.addWidget(lineEdit_write_temp_p, 14, 2, 1, 1)
            frame_26 = QtWidgets.QFrame(frame_2)
            frame_26.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_26.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_26.setObjectName("frame_26")
            gridLayout_write_parameters_entry.addWidget(frame_26, 14, 5, 1, 1)
            label_write_k7 = QtWidgets.QLabel(frame_2)
            label_write_k7.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k7.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k7.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k7.setObjectName("label_write_k7")
            gridLayout_write_parameters_entry.addWidget(label_write_k7, 7, 0, 1, 1)
            frame_28 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_28.sizePolicy().hasHeightForWidth())
            frame_28.setSizePolicy(sizePolicy)
            frame_28.setMinimumSize(QtCore.QSize(10, 0))
            frame_28.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_28.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_28.setObjectName("frame_28")
            gridLayout_write_parameters_entry.addWidget(frame_28, 13, 1, 1, 7)
            label_write_temp = QtWidgets.QLabel(frame_2)
            label_write_temp.setMinimumSize(QtCore.QSize(0, 25))
            label_write_temp.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_temp.setAlignment(QtCore.Qt.AlignCenter)
            label_write_temp.setObjectName("label_write_temp")
            gridLayout_write_parameters_entry.addWidget(label_write_temp, 14, 0, 1, 1)
            frame_23 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_23.sizePolicy().hasHeightForWidth())
            frame_23.setSizePolicy(sizePolicy)
            frame_23.setMinimumSize(QtCore.QSize(10, 0))
            frame_23.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_23.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_23.setObjectName("frame_23")
            gridLayout_write_parameters_entry.addWidget(frame_23, 14, 1, 1, 1)
            label_write_k6 = QtWidgets.QLabel(frame_2)
            label_write_k6.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k6.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k6.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k6.setObjectName("label_write_k6")
            gridLayout_write_parameters_entry.addWidget(label_write_k6, 6, 0, 1, 1)
            label_write_k5 = QtWidgets.QLabel(frame_2)
            label_write_k5.setMinimumSize(QtCore.QSize(0, 25))
            label_write_k5.setMaximumSize(QtCore.QSize(16777215, 16777215))
            label_write_k5.setAlignment(QtCore.Qt.AlignCenter)
            label_write_k5.setObjectName("label_write_k5")
            gridLayout_write_parameters_entry.addWidget(label_write_k5, 5, 0, 1, 1)
            frame_14 = QtWidgets.QFrame(frame_2)
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(frame_14.sizePolicy().hasHeightForWidth())
            frame_14.setSizePolicy(sizePolicy)
            frame_14.setMinimumSize(QtCore.QSize(10, 0))
            frame_14.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame_14.setFrameShadow(QtWidgets.QFrame.Raised)
            frame_14.setObjectName("frame_14")
            gridLayout_write_parameters_entry.addWidget(frame_14, 0, 6, 1, 2)
            verticalLayout.addLayout(gridLayout_write_parameters_entry)
            frame = QtWidgets.QFrame(frame_2)
            frame.setMinimumSize(QtCore.QSize(0, 20))
            frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
            frame.setFrameShadow(QtWidgets.QFrame.Raised)
            frame.setObjectName("frame")
            verticalLayout.addWidget(frame)
            horizontalLayout_save_cancel = QtWidgets.QHBoxLayout()
            horizontalLayout_save_cancel.setObjectName("horizontalLayout_save_cancel")
            buttonBox_write_edit = QtWidgets.QDialogButtonBox(frame_2)
            font = QtGui.QFont()
            font.setKerning(True)
            buttonBox_write_edit.setFont(font)
            buttonBox_write_edit.setOrientation(QtCore.Qt.Horizontal)
            buttonBox_write_edit.setStandardButtons(QtWidgets.QDialogButtonBox.Save)
            buttonBox_write_edit.setObjectName("buttonBox_write_edit")
            horizontalLayout_save_cancel.addWidget(buttonBox_write_edit)
            pushButton_cancel = QtWidgets.QPushButton(frame_2)
            pushButton_cancel.setObjectName("pushButton_cancel")
            horizontalLayout_save_cancel.addWidget(pushButton_cancel)
            verticalLayout.addLayout(horizontalLayout_save_cancel)
            gridLayout.addLayout(verticalLayout, 0, 0, 1, 1)
            gridLayout_2.addWidget(frame_2, 0, 0, 1, 1)

            def retranslateUi(Dialog):
                _translate = QtCore.QCoreApplication.translate
                Dialog.setWindowTitle(_translate("Dialog", "Edit kinetic values"))
                label_write_vdna_defaults.setText(_translate("Dialog", "VDNA defaults"))
                label_write_k4.setText(_translate("Dialog", "k 4"))
                label_write_k8.setText(_translate("Dialog", "k 8"))
                label_write_k3.setText(_translate("Dialog", "k 3"))
                label_write_k9.setText(_translate("Dialog", "k 9"))
                label_write_k11.setText(_translate("Dialog", "k 11"))
                label_write_k2.setText(_translate("Dialog", "k 2"))
                label_write_present_values.setText(_translate("Dialog", "Present values"))
                label_write_tempereture.setText(_translate("Dialog", "Tempereture"))
                pushButton_set_vdna_defaults.setText(_translate("Dialog", "Set VDNA defaults"))
                pushButton_7.setText(_translate("Dialog", "????"))
                label_write_k1.setText(_translate("Dialog", "k 1"))
                pushButton_reset_to_originals.setText(_translate("Dialog", "Reset to originals"))
                label_write_kinetic_values.setText(_translate("Dialog", "Kinetic values"))
                label_write_k10.setText(_translate("Dialog", "k 10"))
                label_write_k7.setText(_translate("Dialog", "k 7"))
                label_write_temp.setText(_translate("Dialog", "Temp"))
                label_write_k6.setText(_translate("Dialog", "k 6"))
                label_write_k5.setText(_translate("Dialog", "k 5"))
                pushButton_cancel.setText(_translate("Dialog", "Cancel"))

            def close_dialog():
                write_advanced_dialog.reject()

            retranslateUi(write_advanced_dialog)
            buttonBox_write_edit.rejected.connect(write_advanced_dialog.reject)
            buttonBox_write_edit.accepted.connect(write_advanced_dialog.accept)

            kt_cells = {'Temp_p': lineEdit_write_temp_p,
                        'k1_p': lineEdit_write_p_k1,
                        'k2_P': lineEdit_write_p_k2,
                        'k3_p': lineEdit_write_p_k3,
                        'k4_p': lineEdit_write_p_k4,
                        'k5_p': lineEdit_write_p_k5,
                        'k6_p': lineEdit_write_p_k6,
                        'k7_p': lineEdit_write_p_k7,
                        'k8_p': lineEdit_write_p_k8,
                        'k9_P': lineEdit_write_p_k9,
                        'k10_p': lineEdit_write_p_k10,
                        'k11_p': lineEdit_write_p_k11}

            set_k_values(source_file, 'get')
            pushButton_set_vdna_defaults.clicked.connect(lambda: set_k_values(source_file, 'set_d'))
            pushButton_reset_to_originals.clicked.connect(lambda: set_k_values(source_file, 'set_o'))
            pushButton_cancel.clicked.connect(close_dialog)

            lineEdit_write_temp_p.textChanged.connect(lambda: validate_cell_values('Temp_p'))
            lineEdit_write_p_k1.textChanged.connect(lambda: validate_cell_values('k1_p'))
            lineEdit_write_p_k2.textChanged.connect(lambda: validate_cell_values('k2_P'))
            lineEdit_write_p_k3.textChanged.connect(lambda: validate_cell_values('k3_p'))
            lineEdit_write_p_k4.textChanged.connect(lambda: validate_cell_values('k4_p'))
            lineEdit_write_p_k5.textChanged.connect(lambda: validate_cell_values('k5_p'))
            lineEdit_write_p_k6.textChanged.connect(lambda: validate_cell_values('k6_p'))
            lineEdit_write_p_k7.textChanged.connect(lambda: validate_cell_values('k7_p'))
            lineEdit_write_p_k8.textChanged.connect(lambda: validate_cell_values('k8_p'))
            lineEdit_write_p_k9.textChanged.connect(lambda: validate_cell_values('k9_P'))
            lineEdit_write_p_k10.textChanged.connect(lambda: validate_cell_values('k10_p'))
            lineEdit_write_p_k11.textChanged.connect(lambda: validate_cell_values('k11_p'))

            run_write_advanced_dialog = write_advanced_dialog.exec_()

            if run_write_advanced_dialog == 1:
                global write_kt_status_g

                write_kt_status_g = True
                save_kt_data()
                self.checkBox_write_advanced.setChecked(True)

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(955, 778)
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(9)
        MainWindow.setFont(font)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.frame_spacing_2 = QtWidgets.QFrame(self.centralwidget)
        self.frame_spacing_2.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_spacing_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_spacing_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_spacing_2.setObjectName("frame_spacing_2")
        self.gridLayout.addWidget(self.frame_spacing_2, 7, 0, 1, 1)
        self.frame_spacing_1 = QtWidgets.QFrame(self.centralwidget)
        self.frame_spacing_1.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_spacing_1.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_spacing_1.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_spacing_1.setObjectName("frame_spacing_1")
        self.gridLayout.addWidget(self.frame_spacing_1, 2, 0, 1, 1)
        self.label_main_title = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(15)
        font.setBold(True)
        font.setWeight(75)
        font.setKerning(True)
        self.label_main_title.setFont(font)
        self.label_main_title.setAlignment(QtCore.Qt.AlignCenter)
        self.label_main_title.setObjectName("label_main_title")
        self.gridLayout.addWidget(self.label_main_title, 0, 0, 1, 1)
        self.label_messages = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        font.setKerning(True)
        self.label_messages.setFont(font)
        self.label_messages.setText("")
        self.label_messages.setAlignment(QtCore.Qt.AlignCenter)
        self.label_messages.setObjectName("label_messages")
        self.gridLayout.addWidget(self.label_messages, 5, 0, 1, 1)
        self.main_tab = QtWidgets.QTabWidget(self.centralwidget)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.main_tab.setFont(font)
        self.main_tab.setTabShape(QtWidgets.QTabWidget.Triangular)
        self.main_tab.setObjectName("main_tab")
        self.tab_write_bngl = QtWidgets.QWidget()
        self.tab_write_bngl.setObjectName("tab_write_bngl")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.tab_write_bngl)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.write_main_frame = QtWidgets.QFrame(self.tab_write_bngl)
        self.write_main_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.write_main_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.write_main_frame.setObjectName("write_main_frame")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.write_main_frame)
        self.verticalLayout.setContentsMargins(50, 50, 50, 50)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout_write_browse_source = QtWidgets.QGridLayout()
        self.gridLayout_write_browse_source.setObjectName("gridLayout_write_browse_source")
        self.label_write_browse_source_info = QtWidgets.QLabel(self.write_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_write_browse_source_info.setFont(font)
        self.label_write_browse_source_info.setAlignment(QtCore.Qt.AlignCenter)
        self.label_write_browse_source_info.setObjectName("label_write_browse_source_info")
        self.gridLayout_write_browse_source.addWidget(self.label_write_browse_source_info, 3, 0, 1, 1)
        self.lineEdit_write_browse_source = QtWidgets.QLineEdit(self.write_main_frame)
        self.lineEdit_write_browse_source.setObjectName("lineEdit_write_browse_source")
        self.gridLayout_write_browse_source.addWidget(self.lineEdit_write_browse_source, 2, 0, 1, 1)
        self.pushButton_write_browse_source = QtWidgets.QPushButton(self.write_main_frame)
        self.pushButton_write_browse_source.setObjectName("pushButton_write_browse_source")
        self.gridLayout_write_browse_source.addWidget(self.pushButton_write_browse_source, 2, 1, 1, 1)
        self.frame_13 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_13.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_13.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_13.setObjectName("frame_13")
        self.gridLayout_write_browse_source.addWidget(self.frame_13, 3, 1, 1, 1)
        self.label_write_browse_source = QtWidgets.QLabel(self.write_main_frame)
        self.label_write_browse_source.setObjectName("label_write_browse_source")
        self.gridLayout_write_browse_source.addWidget(self.label_write_browse_source, 0, 0, 2, 1)
        self.frame_5 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_5.setObjectName("frame_5")
        self.gridLayout_write_browse_source.addWidget(self.frame_5, 0, 1, 2, 1)
        self.verticalLayout.addLayout(self.gridLayout_write_browse_source)
        self.frame_10 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_10.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_10.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_10.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_10.setObjectName("frame_10")
        self.verticalLayout.addWidget(self.frame_10)
        self.verticalLayout_write_list_created = QtWidgets.QVBoxLayout()
        self.verticalLayout_write_list_created.setObjectName("verticalLayout_write_list_created")
        self.label_write_list_created = QtWidgets.QLabel(self.write_main_frame)
        self.label_write_list_created.setObjectName("label_write_list_created")
        self.verticalLayout_write_list_created.addWidget(self.label_write_list_created)
        self.listWidget_write_list_created = QtWidgets.QListWidget(self.write_main_frame)
        self.listWidget_write_list_created.setMinimumSize(QtCore.QSize(0, 0))
        font = QtGui.QFont()
        font.setFamily("Lucida Console")
        self.listWidget_write_list_created.setFont(font)
        self.listWidget_write_list_created.setObjectName("listWidget_write_list_created")
        self.verticalLayout_write_list_created.addWidget(self.listWidget_write_list_created)
        self.verticalLayout.addLayout(self.verticalLayout_write_list_created)
        self.horizontalLayout_write_tool_buttons = QtWidgets.QHBoxLayout()
        self.horizontalLayout_write_tool_buttons.setObjectName("horizontalLayout_write_tool_buttons")
        self.pushButton_write_edit = QtWidgets.QPushButton(self.write_main_frame)
        self.pushButton_write_edit.setObjectName("pushButton_write_edit")
        self.horizontalLayout_write_tool_buttons.addWidget(self.pushButton_write_edit)
        self.frame_19 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_19.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_19.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_19.setObjectName("frame_19")
        self.horizontalLayout_write_tool_buttons.addWidget(self.frame_19)
        self.pushButton_write_delete = QtWidgets.QPushButton(self.write_main_frame)
        self.pushButton_write_delete.setObjectName("pushButton_write_delete")
        self.horizontalLayout_write_tool_buttons.addWidget(self.pushButton_write_delete)
        self.frame_20 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_20.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_20.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_20.setObjectName("frame_20")
        self.horizontalLayout_write_tool_buttons.addWidget(self.frame_20)
        self.pushButton_write_reset_all = QtWidgets.QPushButton(self.write_main_frame)
        self.pushButton_write_reset_all.setObjectName("pushButton_write_reset_all")
        self.horizontalLayout_write_tool_buttons.addWidget(self.pushButton_write_reset_all)
        self.verticalLayout.addLayout(self.horizontalLayout_write_tool_buttons)
        self.frame_35 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_35.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_35.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_35.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_35.setObjectName("frame_35")
        self.verticalLayout.addWidget(self.frame_35)
        self.gridLayout_write_sequence_entry = QtWidgets.QGridLayout()
        self.gridLayout_write_sequence_entry.setObjectName("gridLayout_write_sequence_entry")
        self.frame_7 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_7.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_7.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_7.setObjectName("frame_7")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_7, 0, 1, 1, 1)
        self.lineEdit_custom_file_name = QtWidgets.QLineEdit(self.write_main_frame)
        self.lineEdit_custom_file_name.setEnabled(False)
        self.lineEdit_custom_file_name.setObjectName("lineEdit_custom_file_name")
        self.gridLayout_write_sequence_entry.addWidget(self.lineEdit_custom_file_name, 4, 1, 1, 3)
        self.checkBox_custom_file_name = QtWidgets.QCheckBox(self.write_main_frame)
        self.checkBox_custom_file_name.setObjectName("checkBox_custom_file_name")
        self.gridLayout_write_sequence_entry.addWidget(self.checkBox_custom_file_name, 4, 0, 1, 1)
        self.pushButton_write_import_sequence = QtWidgets.QPushButton(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_write_import_sequence.sizePolicy().hasHeightForWidth())
        self.pushButton_write_import_sequence.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_write_import_sequence.setFont(font)
        self.pushButton_write_import_sequence.setObjectName("pushButton_write_import_sequence")
        self.gridLayout_write_sequence_entry.addWidget(self.pushButton_write_import_sequence, 2, 4, 1, 1)
        self.frame_43 = QtWidgets.QFrame(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_43.sizePolicy().hasHeightForWidth())
        self.frame_43.setSizePolicy(sizePolicy)
        self.frame_43.setMinimumSize(QtCore.QSize(0, 5))
        self.frame_43.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_43.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_43.setObjectName("frame_43")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_43, 5, 0, 1, 5)
        self.pushButton_write_submit_without = QtWidgets.QPushButton(self.write_main_frame)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_write_submit_without.setFont(font)
        self.pushButton_write_submit_without.setObjectName("pushButton_write_submit_without")
        self.gridLayout_write_sequence_entry.addWidget(self.pushButton_write_submit_without, 6, 3, 1, 2)
        self.frame_6 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_6.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_6.setObjectName("frame_6")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_6, 6, 2, 1, 1)
        self.label_write_allowed_char = QtWidgets.QLabel(self.write_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_write_allowed_char.setFont(font)
        self.label_write_allowed_char.setAlignment(QtCore.Qt.AlignCenter)
        self.label_write_allowed_char.setObjectName("label_write_allowed_char")
        self.gridLayout_write_sequence_entry.addWidget(self.label_write_allowed_char, 2, 1, 1, 1)
        self.pushButton_write_submit_with = QtWidgets.QPushButton(self.write_main_frame)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_write_submit_with.setFont(font)
        self.pushButton_write_submit_with.setObjectName("pushButton_write_submit_with")
        self.gridLayout_write_sequence_entry.addWidget(self.pushButton_write_submit_with, 6, 0, 1, 2)
        self.frame_30 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_30.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_30.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_30.setObjectName("frame_30")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_30, 2, 2, 1, 1)
        self.label_write_amount = QtWidgets.QLabel(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_write_amount.sizePolicy().hasHeightForWidth())
        self.label_write_amount.setSizePolicy(sizePolicy)
        self.label_write_amount.setAlignment(QtCore.Qt.AlignCenter)
        self.label_write_amount.setObjectName("label_write_amount")
        self.gridLayout_write_sequence_entry.addWidget(self.label_write_amount, 0, 3, 1, 1)
        self.frame_32 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_32.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_32.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_32.setObjectName("frame_32")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_32, 2, 0, 1, 1)
        self.frame_4 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4.setObjectName("frame_4")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_4, 0, 4, 1, 1)
        self.lineEdit_write_sequence = QtWidgets.QLineEdit(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_write_sequence.sizePolicy().hasHeightForWidth())
        self.lineEdit_write_sequence.setSizePolicy(sizePolicy)
        self.lineEdit_write_sequence.setObjectName("lineEdit_write_sequence")
        self.gridLayout_write_sequence_entry.addWidget(self.lineEdit_write_sequence, 1, 0, 1, 3)
        self.frame_8 = QtWidgets.QFrame(self.write_main_frame)
        self.frame_8.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_8.setObjectName("frame_8")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_8, 0, 2, 1, 1)
        self.lineEdit_write_amount = QtWidgets.QLineEdit(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_write_amount.sizePolicy().hasHeightForWidth())
        self.lineEdit_write_amount.setSizePolicy(sizePolicy)
        self.lineEdit_write_amount.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_write_amount.setObjectName("lineEdit_write_amount")
        self.gridLayout_write_sequence_entry.addWidget(self.lineEdit_write_amount, 1, 3, 1, 1)
        self.label_write_allowed_num = QtWidgets.QLabel(self.write_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_write_allowed_num.setFont(font)
        self.label_write_allowed_num.setAlignment(QtCore.Qt.AlignCenter)
        self.label_write_allowed_num.setObjectName("label_write_allowed_num")
        self.gridLayout_write_sequence_entry.addWidget(self.label_write_allowed_num, 2, 3, 1, 1)
        self.pushButton_write_add_sequence = QtWidgets.QPushButton(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_write_add_sequence.sizePolicy().hasHeightForWidth())
        self.pushButton_write_add_sequence.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_write_add_sequence.setFont(font)
        self.pushButton_write_add_sequence.setObjectName("pushButton_write_add_sequence")
        self.gridLayout_write_sequence_entry.addWidget(self.pushButton_write_add_sequence, 1, 4, 1, 1)
        self.label_write_sequence = QtWidgets.QLabel(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_write_sequence.sizePolicy().hasHeightForWidth())
        self.label_write_sequence.setSizePolicy(sizePolicy)
        self.label_write_sequence.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_write_sequence.setObjectName("label_write_sequence")
        self.gridLayout_write_sequence_entry.addWidget(self.label_write_sequence, 0, 0, 1, 1)
        self.frame_34 = QtWidgets.QFrame(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_34.sizePolicy().hasHeightForWidth())
        self.frame_34.setSizePolicy(sizePolicy)
        self.frame_34.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_34.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_34.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_34.setObjectName("frame_34")
        self.gridLayout_write_sequence_entry.addWidget(self.frame_34, 3, 0, 1, 5)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.checkBox_write_advanced = QtWidgets.QCheckBox(self.write_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_write_advanced.sizePolicy().hasHeightForWidth())
        self.checkBox_write_advanced.setSizePolicy(sizePolicy)
        self.checkBox_write_advanced.setText("")
        self.checkBox_write_advanced.setObjectName("checkBox_write_advanced")
        self.horizontalLayout_2.addWidget(self.checkBox_write_advanced)
        self.pushButton_write_advanced = QtWidgets.QPushButton(self.write_main_frame)
        self.pushButton_write_advanced.setObjectName("pushButton_write_advanced")
        self.horizontalLayout_2.addWidget(self.pushButton_write_advanced)
        self.gridLayout_write_sequence_entry.addLayout(self.horizontalLayout_2, 4, 4, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout_write_sequence_entry)
        self.gridLayout_2.addWidget(self.write_main_frame, 0, 0, 1, 1)
        self.main_tab.addTab(self.tab_write_bngl, "")
        self.tab_run_bngl = QtWidgets.QWidget()
        self.tab_run_bngl.setObjectName("tab_run_bngl")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.tab_run_bngl)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.run_main_frame = QtWidgets.QFrame(self.tab_run_bngl)
        self.run_main_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.run_main_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.run_main_frame.setObjectName("run_main_frame")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.run_main_frame)
        self.verticalLayout_3.setContentsMargins(50, 50, 50, 50)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.gridLayout_run_browse_source = QtWidgets.QGridLayout()
        self.gridLayout_run_browse_source.setObjectName("gridLayout_run_browse_source")
        self.label_run_browse_source_info = QtWidgets.QLabel(self.run_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_run_browse_source_info.setFont(font)
        self.label_run_browse_source_info.setAlignment(QtCore.Qt.AlignCenter)
        self.label_run_browse_source_info.setObjectName("label_run_browse_source_info")
        self.gridLayout_run_browse_source.addWidget(self.label_run_browse_source_info, 3, 0, 1, 1)
        self.pushButton_run_browse_source = QtWidgets.QPushButton(self.run_main_frame)
        self.pushButton_run_browse_source.setObjectName("pushButton_run_browse_source")
        self.gridLayout_run_browse_source.addWidget(self.pushButton_run_browse_source, 2, 1, 1, 1)
        self.frame_16 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_16.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_16.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_16.setObjectName("frame_16")
        self.gridLayout_run_browse_source.addWidget(self.frame_16, 3, 1, 1, 1)
        self.lineEdit_run_browse_source = QtWidgets.QLineEdit(self.run_main_frame)
        self.lineEdit_run_browse_source.setObjectName("lineEdit_run_browse_source")
        self.gridLayout_run_browse_source.addWidget(self.lineEdit_run_browse_source, 2, 0, 1, 1)
        self.label_run_browse_source = QtWidgets.QLabel(self.run_main_frame)
        self.label_run_browse_source.setObjectName("label_run_browse_source")
        self.gridLayout_run_browse_source.addWidget(self.label_run_browse_source, 0, 0, 2, 1)
        self.frame_15 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_15.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_15.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_15.setObjectName("frame_15")
        self.gridLayout_run_browse_source.addWidget(self.frame_15, 0, 1, 2, 1)
        self.verticalLayout_3.addLayout(self.gridLayout_run_browse_source)
        self.frame = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame.sizePolicy().hasHeightForWidth())
        self.frame.setSizePolicy(sizePolicy)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout_3.addWidget(self.frame)
        self.gridLayout_run_browse_source_des = QtWidgets.QGridLayout()
        self.gridLayout_run_browse_source_des.setObjectName("gridLayout_run_browse_source_des")
        self.pushButton_run_browse_des = QtWidgets.QPushButton(self.run_main_frame)
        self.pushButton_run_browse_des.setObjectName("pushButton_run_browse_des")
        self.gridLayout_run_browse_source_des.addWidget(self.pushButton_run_browse_des, 2, 1, 1, 1)
        self.label_run_browse_des_info = QtWidgets.QLabel(self.run_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_run_browse_des_info.setFont(font)
        self.label_run_browse_des_info.setAlignment(QtCore.Qt.AlignCenter)
        self.label_run_browse_des_info.setObjectName("label_run_browse_des_info")
        self.gridLayout_run_browse_source_des.addWidget(self.label_run_browse_des_info, 3, 0, 1, 1)
        self.lineEdit_run_browse_des = QtWidgets.QLineEdit(self.run_main_frame)
        self.lineEdit_run_browse_des.setObjectName("lineEdit_run_browse_des")
        self.gridLayout_run_browse_source_des.addWidget(self.lineEdit_run_browse_des, 2, 0, 1, 1)
        self.frame_18 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_18.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_18.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_18.setObjectName("frame_18")
        self.gridLayout_run_browse_source_des.addWidget(self.frame_18, 3, 1, 1, 1)
        self.label_run_browse_des = QtWidgets.QLabel(self.run_main_frame)
        self.label_run_browse_des.setObjectName("label_run_browse_des")
        self.gridLayout_run_browse_source_des.addWidget(self.label_run_browse_des, 0, 0, 2, 1)
        self.frame_17 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_17.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_17.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_17.setObjectName("frame_17")
        self.gridLayout_run_browse_source_des.addWidget(self.frame_17, 0, 1, 2, 1)
        self.verticalLayout_3.addLayout(self.gridLayout_run_browse_source_des)
        self.frame_2 = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.verticalLayout_3.addWidget(self.frame_2)
        self.gridLayout_run_parameters_entry = QtWidgets.QGridLayout()
        self.gridLayout_run_parameters_entry.setObjectName("gridLayout_run_parameters_entry")
        self.frame_47 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_47.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_47.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_47.setObjectName("frame_47")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_47, 0, 6, 1, 1)
        self.pushButton_run_advanced = QtWidgets.QPushButton(self.run_main_frame)
        self.pushButton_run_advanced.setObjectName("pushButton_run_advanced")
        self.gridLayout_run_parameters_entry.addWidget(self.pushButton_run_advanced, 1, 6, 1, 1)
        self.checkBox_run_advanced = QtWidgets.QCheckBox(self.run_main_frame)
        self.checkBox_run_advanced.setText("")
        self.checkBox_run_advanced.setObjectName("checkBox_run_advanced")
        self.gridLayout_run_parameters_entry.addWidget(self.checkBox_run_advanced, 1, 7, 1, 1)
        self.frame_48 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_48.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_48.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_48.setObjectName("frame_48")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_48, 0, 7, 1, 1)
        self.lineEdit_run_n_dumps = QtWidgets.QLineEdit(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_run_n_dumps.sizePolicy().hasHeightForWidth())
        self.lineEdit_run_n_dumps.setSizePolicy(sizePolicy)
        self.lineEdit_run_n_dumps.setMinimumSize(QtCore.QSize(100, 0))
        self.lineEdit_run_n_dumps.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_run_n_dumps.setObjectName("lineEdit_run_n_dumps")
        self.gridLayout_run_parameters_entry.addWidget(self.lineEdit_run_n_dumps, 1, 2, 1, 1)
        self.frame_12 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_12.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_12.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_12.setObjectName("frame_12")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_12, 0, 9, 1, 2)
        self.frame_11 = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_11.sizePolicy().hasHeightForWidth())
        self.frame_11.setSizePolicy(sizePolicy)
        self.frame_11.setMinimumSize(QtCore.QSize(10, 0))
        self.frame_11.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_11.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_11.setObjectName("frame_11")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_11, 0, 3, 2, 1)
        self.lineEdit_run_start_time = QtWidgets.QLineEdit(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_run_start_time.sizePolicy().hasHeightForWidth())
        self.lineEdit_run_start_time.setSizePolicy(sizePolicy)
        self.lineEdit_run_start_time.setMinimumSize(QtCore.QSize(100, 0))
        self.lineEdit_run_start_time.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_run_start_time.setObjectName("lineEdit_run_start_time")
        self.gridLayout_run_parameters_entry.addWidget(self.lineEdit_run_start_time, 1, 0, 1, 1)
        self.frame_24 = QtWidgets.QFrame(self.run_main_frame)
        self.frame_24.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_24.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_24.setObjectName("frame_24")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_24, 2, 0, 1, 11)
        self.label_run_start_time = QtWidgets.QLabel(self.run_main_frame)
        self.label_run_start_time.setAlignment(QtCore.Qt.AlignCenter)
        self.label_run_start_time.setObjectName("label_run_start_time")
        self.gridLayout_run_parameters_entry.addWidget(self.label_run_start_time, 0, 0, 1, 1)
        self.frame_25 = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_25.sizePolicy().hasHeightForWidth())
        self.frame_25.setSizePolicy(sizePolicy)
        self.frame_25.setMinimumSize(QtCore.QSize(50, 0))
        self.frame_25.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_25.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_25.setObjectName("frame_25")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_25, 0, 5, 2, 1)
        self.label_run_n_dumps = QtWidgets.QLabel(self.run_main_frame)
        self.label_run_n_dumps.setAlignment(QtCore.Qt.AlignCenter)
        self.label_run_n_dumps.setObjectName("label_run_n_dumps")
        self.gridLayout_run_parameters_entry.addWidget(self.label_run_n_dumps, 0, 2, 1, 1)
        self.label_run_sim_end = QtWidgets.QLabel(self.run_main_frame)
        self.label_run_sim_end.setAlignment(QtCore.Qt.AlignCenter)
        self.label_run_sim_end.setObjectName("label_run_sim_end")
        self.gridLayout_run_parameters_entry.addWidget(self.label_run_sim_end, 0, 4, 1, 1)
        self.frame_9 = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_9.sizePolicy().hasHeightForWidth())
        self.frame_9.setSizePolicy(sizePolicy)
        self.frame_9.setMinimumSize(QtCore.QSize(10, 0))
        self.frame_9.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_9.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_9.setObjectName("frame_9")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_9, 0, 1, 2, 1)
        self.pushButton_run_run = QtWidgets.QPushButton(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_run_run.sizePolicy().hasHeightForWidth())
        self.pushButton_run_run.setSizePolicy(sizePolicy)
        self.pushButton_run_run.setMinimumSize(QtCore.QSize(150, 0))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_run_run.setFont(font)
        self.pushButton_run_run.setObjectName("pushButton_run_run")
        self.gridLayout_run_parameters_entry.addWidget(self.pushButton_run_run, 1, 9, 1, 2)
        self.lineEdit_run_sim_end = QtWidgets.QLineEdit(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_run_sim_end.sizePolicy().hasHeightForWidth())
        self.lineEdit_run_sim_end.setSizePolicy(sizePolicy)
        self.lineEdit_run_sim_end.setMinimumSize(QtCore.QSize(100, 0))
        self.lineEdit_run_sim_end.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit_run_sim_end.setObjectName("lineEdit_run_sim_end")
        self.gridLayout_run_parameters_entry.addWidget(self.lineEdit_run_sim_end, 1, 4, 1, 1)
        self.frame_46 = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_46.sizePolicy().hasHeightForWidth())
        self.frame_46.setSizePolicy(sizePolicy)
        self.frame_46.setMinimumSize(QtCore.QSize(50, 0))
        self.frame_46.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_46.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_46.setObjectName("frame_46")
        self.gridLayout_run_parameters_entry.addWidget(self.frame_46, 0, 8, 2, 1)
        self.verticalLayout_3.addLayout(self.gridLayout_run_parameters_entry)
        self.frame_3 = QtWidgets.QFrame(self.run_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_3.sizePolicy().hasHeightForWidth())
        self.frame_3.setSizePolicy(sizePolicy)
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.verticalLayout_3.addWidget(self.frame_3)
        self.gridLayout_4.addWidget(self.run_main_frame, 1, 0, 1, 1)
        self.main_tab.addTab(self.tab_run_bngl, "")
        self.tab_visualize = QtWidgets.QWidget()
        self.tab_visualize.setObjectName("tab_visualize")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.tab_visualize)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.visualize_main_frame = QtWidgets.QFrame(self.tab_visualize)
        self.visualize_main_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.visualize_main_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.visualize_main_frame.setObjectName("visualize_main_frame")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.visualize_main_frame)
        self.verticalLayout_2.setContentsMargins(50, 50, 50, 50)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gridLayout_visual_browse_source = QtWidgets.QGridLayout()
        self.gridLayout_visual_browse_source.setObjectName("gridLayout_visual_browse_source")
        self.label_visual_browse_source_info = QtWidgets.QLabel(self.visualize_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_visual_browse_source_info.setFont(font)
        self.label_visual_browse_source_info.setAlignment(QtCore.Qt.AlignCenter)
        self.label_visual_browse_source_info.setObjectName("label_visual_browse_source_info")
        self.gridLayout_visual_browse_source.addWidget(self.label_visual_browse_source_info, 3, 0, 1, 1)
        self.lineEdit_visual_browse_source = QtWidgets.QLineEdit(self.visualize_main_frame)
        self.lineEdit_visual_browse_source.setObjectName("lineEdit_visual_browse_source")
        self.gridLayout_visual_browse_source.addWidget(self.lineEdit_visual_browse_source, 2, 0, 1, 1)
        self.pushButton_visual_browse_source = QtWidgets.QPushButton(self.visualize_main_frame)
        self.pushButton_visual_browse_source.setObjectName("pushButton_visual_browse_source")
        self.gridLayout_visual_browse_source.addWidget(self.pushButton_visual_browse_source, 2, 1, 1, 1)
        self.frame_23 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_23.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_23.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_23.setObjectName("frame_23")
        self.gridLayout_visual_browse_source.addWidget(self.frame_23, 3, 1, 1, 1)
        self.label_visual_browse_source = QtWidgets.QLabel(self.visualize_main_frame)
        self.label_visual_browse_source.setObjectName("label_visual_browse_source")
        self.gridLayout_visual_browse_source.addWidget(self.label_visual_browse_source, 0, 0, 2, 1)
        self.frame_27 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_27.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_27.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_27.setObjectName("frame_27")
        self.gridLayout_visual_browse_source.addWidget(self.frame_27, 0, 1, 2, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_visual_browse_source)
        self.frame_28 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_28.sizePolicy().hasHeightForWidth())
        self.frame_28.setSizePolicy(sizePolicy)
        self.frame_28.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_28.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_28.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_28.setObjectName("frame_28")
        self.verticalLayout_2.addWidget(self.frame_28)
        self.gridLayout_visual_history = QtWidgets.QGridLayout()
        self.gridLayout_visual_history.setObjectName("gridLayout_visual_history")
        self.pushButton_visual_view_source = QtWidgets.QPushButton(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_visual_view_source.sizePolicy().hasHeightForWidth())
        self.pushButton_visual_view_source.setSizePolicy(sizePolicy)
        self.pushButton_visual_view_source.setMinimumSize(QtCore.QSize(200, 0))
        self.pushButton_visual_view_source.setObjectName("pushButton_visual_view_source")
        self.gridLayout_visual_history.addWidget(self.pushButton_visual_view_source, 6, 2, 1, 1)
        self.frame_39 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_39.sizePolicy().hasHeightForWidth())
        self.frame_39.setSizePolicy(sizePolicy)
        self.frame_39.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_39.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_39.setObjectName("frame_39")
        self.gridLayout_visual_history.addWidget(self.frame_39, 7, 2, 1, 1)
        self.pushButton_visual_delete = QtWidgets.QPushButton(self.visualize_main_frame)
        self.pushButton_visual_delete.setMinimumSize(QtCore.QSize(200, 0))
        self.pushButton_visual_delete.setObjectName("pushButton_visual_delete")
        self.gridLayout_visual_history.addWidget(self.pushButton_visual_delete, 8, 2, 1, 1)
        self.pushButton_visual_view_bngl = QtWidgets.QPushButton(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_visual_view_bngl.sizePolicy().hasHeightForWidth())
        self.pushButton_visual_view_bngl.setSizePolicy(sizePolicy)
        self.pushButton_visual_view_bngl.setMinimumSize(QtCore.QSize(200, 0))
        self.pushButton_visual_view_bngl.setObjectName("pushButton_visual_view_bngl")
        self.gridLayout_visual_history.addWidget(self.pushButton_visual_view_bngl, 4, 2, 1, 1)
        self.pushButton_visual_view_comp = QtWidgets.QPushButton(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_visual_view_comp.sizePolicy().hasHeightForWidth())
        self.pushButton_visual_view_comp.setSizePolicy(sizePolicy)
        self.pushButton_visual_view_comp.setMinimumSize(QtCore.QSize(200, 0))
        self.pushButton_visual_view_comp.setObjectName("pushButton_visual_view_comp")
        self.gridLayout_visual_history.addWidget(self.pushButton_visual_view_comp, 2, 2, 1, 1)
        self.label_visual_history = QtWidgets.QLabel(self.visualize_main_frame)
        self.label_visual_history.setObjectName("label_visual_history")
        self.gridLayout_visual_history.addWidget(self.label_visual_history, 0, 0, 1, 1)
        self.frame_41 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_41.sizePolicy().hasHeightForWidth())
        self.frame_41.setSizePolicy(sizePolicy)
        self.frame_41.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_41.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_41.setObjectName("frame_41")
        self.gridLayout_visual_history.addWidget(self.frame_41, 9, 2, 1, 1)
        self.frame_42 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_42.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_42.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_42.setObjectName("frame_42")
        self.gridLayout_visual_history.addWidget(self.frame_42, 0, 2, 1, 1)
        self.frame_38 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_38.sizePolicy().hasHeightForWidth())
        self.frame_38.setSizePolicy(sizePolicy)
        self.frame_38.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_38.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_38.setObjectName("frame_38")
        self.gridLayout_visual_history.addWidget(self.frame_38, 5, 2, 1, 1)
        self.frame_29 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_29.sizePolicy().hasHeightForWidth())
        self.frame_29.setSizePolicy(sizePolicy)
        self.frame_29.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_29.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_29.setObjectName("frame_29")
        self.gridLayout_visual_history.addWidget(self.frame_29, 1, 2, 1, 1)
        self.frame_36 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_36.sizePolicy().hasHeightForWidth())
        self.frame_36.setSizePolicy(sizePolicy)
        self.frame_36.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_36.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_36.setObjectName("frame_36")
        self.gridLayout_visual_history.addWidget(self.frame_36, 3, 2, 1, 1)
        self.listWidget_visual_history_list = QtWidgets.QListWidget(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.listWidget_visual_history_list.sizePolicy().hasHeightForWidth())
        self.listWidget_visual_history_list.setSizePolicy(sizePolicy)
        self.listWidget_visual_history_list.setMinimumSize(QtCore.QSize(0, 0))
        self.listWidget_visual_history_list.setObjectName("listWidget_visual_history_list")
        self.gridLayout_visual_history.addWidget(self.listWidget_visual_history_list, 1, 0, 9, 1)
        self.frame_21 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_21.setMinimumSize(QtCore.QSize(50, 0))
        self.frame_21.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_21.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_21.setObjectName("frame_21")
        self.gridLayout_visual_history.addWidget(self.frame_21, 0, 1, 10, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_visual_history)
        self.frame_33 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_33.sizePolicy().hasHeightForWidth())
        self.frame_33.setSizePolicy(sizePolicy)
        self.frame_33.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_33.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_33.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_33.setObjectName("frame_33")
        self.verticalLayout_2.addWidget(self.frame_33)
        self.gridLayout_visual_browse_des = QtWidgets.QGridLayout()
        self.gridLayout_visual_browse_des.setObjectName("gridLayout_visual_browse_des")
        self.pushButton_visual_browse_des = QtWidgets.QPushButton(self.visualize_main_frame)
        self.pushButton_visual_browse_des.setObjectName("pushButton_visual_browse_des")
        self.gridLayout_visual_browse_des.addWidget(self.pushButton_visual_browse_des, 2, 1, 1, 1)
        self.label_visual_browse_des_info = QtWidgets.QLabel(self.visualize_main_frame)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label_visual_browse_des_info.setFont(font)
        self.label_visual_browse_des_info.setAlignment(QtCore.Qt.AlignCenter)
        self.label_visual_browse_des_info.setObjectName("label_visual_browse_des_info")
        self.gridLayout_visual_browse_des.addWidget(self.label_visual_browse_des_info, 3, 0, 1, 1)
        self.lineEdit_visual_browse_des = QtWidgets.QLineEdit(self.visualize_main_frame)
        self.lineEdit_visual_browse_des.setObjectName("lineEdit_visual_browse_des")
        self.gridLayout_visual_browse_des.addWidget(self.lineEdit_visual_browse_des, 2, 0, 1, 1)
        self.frame_44 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_44.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_44.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_44.setObjectName("frame_44")
        self.gridLayout_visual_browse_des.addWidget(self.frame_44, 3, 1, 1, 1)
        self.label_visual_browse_des = QtWidgets.QLabel(self.visualize_main_frame)
        self.label_visual_browse_des.setObjectName("label_visual_browse_des")
        self.gridLayout_visual_browse_des.addWidget(self.label_visual_browse_des, 0, 0, 2, 1)
        self.frame_45 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_45.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_45.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_45.setObjectName("frame_45")
        self.gridLayout_visual_browse_des.addWidget(self.frame_45, 0, 1, 2, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_visual_browse_des)
        self.frame_37 = QtWidgets.QFrame(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_37.sizePolicy().hasHeightForWidth())
        self.frame_37.setSizePolicy(sizePolicy)
        self.frame_37.setMinimumSize(QtCore.QSize(0, 15))
        self.frame_37.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_37.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_37.setObjectName("frame_37")
        self.verticalLayout_2.addWidget(self.frame_37)
        self.gridLayout_visual_run = QtWidgets.QGridLayout()
        self.gridLayout_visual_run.setObjectName("gridLayout_visual_run")
        self.pushButton_visual_advanced = QtWidgets.QPushButton(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_visual_advanced.sizePolicy().hasHeightForWidth())
        self.pushButton_visual_advanced.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_visual_advanced.setFont(font)
        self.pushButton_visual_advanced.setObjectName("pushButton_visual_advanced")
        self.gridLayout_visual_run.addWidget(self.pushButton_visual_advanced, 0, 2, 1, 1)
        self.checkBox_visual_advanced = QtWidgets.QCheckBox(self.visualize_main_frame)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_visual_advanced.sizePolicy().hasHeightForWidth())
        self.checkBox_visual_advanced.setSizePolicy(sizePolicy)
        self.checkBox_visual_advanced.setText("")
        self.checkBox_visual_advanced.setObjectName("checkBox_visual_advanced")
        self.gridLayout_visual_run.addWidget(self.checkBox_visual_advanced, 0, 3, 1, 1)
        self.frame_40 = QtWidgets.QFrame(self.visualize_main_frame)
        self.frame_40.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_40.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_40.setObjectName("frame_40")
        self.gridLayout_visual_run.addWidget(self.frame_40, 0, 4, 1, 2)
        self.pushButton_visual_visualize = QtWidgets.QPushButton(self.visualize_main_frame)
        font = QtGui.QFont()
        font.setPointSize(11)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_visual_visualize.setFont(font)
        self.pushButton_visual_visualize.setObjectName("pushButton_visual_visualize")
        self.gridLayout_visual_run.addWidget(self.pushButton_visual_visualize, 0, 7, 1, 1)
        self.verticalLayout_2.addLayout(self.gridLayout_visual_run)
        self.gridLayout_6.addWidget(self.visualize_main_frame, 0, 0, 1, 1)
        self.main_tab.addTab(self.tab_visualize, "")
        self.gridLayout.addWidget(self.main_tab, 9, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.main_tab.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)


        """ all commands """
        # write bngl tab
        self.lineEdit_write_amount.setText('1')
        self.pushButton_write_browse_source.clicked.connect(lambda: browse_path(write_b_source, write_tab))
        self.pushButton_write_submit_with.clicked.connect(lambda: write_bngl('with'))
        self.pushButton_write_submit_without.clicked.connect(lambda: write_bngl('without'))
        self.lineEdit_write_sequence.textChanged.connect(lambda: validate_char(write_the_seq))
        self.lineEdit_write_amount.textChanged.connect(lambda: validate_num(write_n_seq, write_tab))
        self.lineEdit_write_browse_source.textChanged.connect(lambda: on_edit_source_path(write_b_source, write_tab))
        self.pushButton_write_add_sequence.clicked.connect(add_sequence)
        self.pushButton_write_delete.clicked.connect(
            lambda: delete_created_sequences(self.listWidget_write_list_created.currentRow()))
        self.pushButton_write_reset_all.clicked.connect(reset_created_sequences_list)
        self.listWidget_write_list_created.itemClicked.connect(lambda: validate_and_update_buttons(write_tab))
        self.pushButton_write_reset_all.setDisabled(True)
        self.pushButton_write_delete.setDisabled(True)
        self.pushButton_write_edit.setDisabled(True)
        self.pushButton_write_add_sequence.setDisabled(True)
        self.pushButton_write_submit_with.setDisabled(True)
        self.pushButton_write_submit_without.setDisabled(True)
        self.pushButton_write_import_sequence.clicked.connect(lambda: import_sequences(import_seq, write_tab))
        self.pushButton_write_edit.clicked.connect(lambda: edit_seq(self.listWidget_write_list_created.currentRow()))
        self.checkBox_custom_file_name.clicked.connect(custom_name_set)
        self.lineEdit_custom_file_name.textChanged.connect(validate_custom_name)
        self.pushButton_write_advanced.setDisabled(True)
        self.pushButton_write_advanced.clicked.connect(write_advanced)
        self.checkBox_write_advanced.clicked.connect(write_activate_write_advanced)

        # run bngl
        self.lineEdit_run_start_time.setText('0')
        self.lineEdit_run_n_dumps.setText('10')
        self.lineEdit_run_sim_end.setText('1')
        self.pushButton_run_browse_source.clicked.connect(lambda: browse_path(run_b_source, run_tab))
        self.pushButton_run_browse_des.clicked.connect(lambda: browse_path(run_b_des, run_tab))
        self.pushButton_run_run.setDisabled(True)
        self.lineEdit_run_start_time.textChanged.connect(lambda: validate_num(run_start_time, run_tab))
        self.lineEdit_run_n_dumps.textChanged.connect(lambda: validate_num(run_n_dumps, run_tab))
        self.lineEdit_run_sim_end.textChanged.connect(lambda: validate_num(run_sim_end, run_tab))
        self.lineEdit_run_browse_source.textChanged.connect(lambda: on_edit_source_path(run_b_source, run_tab))
        self.lineEdit_run_browse_des.textChanged.connect(lambda: validate_and_update_buttons(run_tab))
        self.pushButton_run_run.clicked.connect(run_bngl)
        self.pushButton_run_advanced.clicked.connect(run_advanced)
        self.pushButton_run_advanced.setDisabled(True)
        self.checkBox_run_advanced.clicked.connect(run_activate_run_advanced)

        # visual bngl
        self.pushButton_visual_visualize.setDisabled(True)
        self.pushButton_visual_view_comp.setDisabled(True)
        self.pushButton_visual_view_source.setDisabled(True)
        self.pushButton_visual_view_bngl.setDisabled(True)
        self.pushButton_visual_delete.setDisabled(True)
        self.pushButton_visual_visualize.clicked.connect(lambda: run_visualisation(''))
        self.pushButton_visual_view_comp.clicked.connect(lambda: pop_view('.html', 'html'))
        self.pushButton_visual_view_source.clicked.connect(lambda: pop_view('.species', 's_sp'))
        self.pushButton_visual_view_bngl.clicked.connect(lambda: pop_view('.species', 'r_sp'))
        self.pushButton_visual_delete.clicked.connect(delete_folder)
        self.pushButton_visual_browse_source.clicked.connect(lambda: browse_path(visual_b_source, visual_tab))
        self.pushButton_visual_browse_des.clicked.connect(lambda: browse_path(visual_b_des, visual_tab))
        self.lineEdit_visual_browse_source.textChanged.connect(lambda: on_edit_source_path(visual_b_source, visual_tab))
        self.lineEdit_visual_browse_des.textChanged.connect(lambda: validate_and_update_buttons(visual_tab))
        self.listWidget_visual_history_list.itemClicked.connect(lambda: validate_and_update_buttons(visual_tab))
        self.lineEdit_visual_browse_des.textChanged.connect(update_history)
        self.pushButton_visual_advanced.clicked.connect(advanced_options)

        def restore_previous_session():
            all_source_cells = [[write_b_source, self.lineEdit_write_browse_source],
                                [run_b_source, self.lineEdit_run_browse_source],
                                [visual_b_source, self.lineEdit_visual_browse_source]]

            all_des_cells = [[run_b_des, self.lineEdit_run_browse_des],
                             [visual_b_des, self.lineEdit_visual_browse_des]]

            for s_cells in all_source_cells:
                file_link = read_s_loc(s_cells[0])
                if os.path.isfile(file_link):
                    s_cells[1].setText(file_link)

            for des_cells in all_des_cells:
                dir_link = read_s_loc(des_cells[0])
                if os.path.isdir(dir_link):
                    des_cells[1].setText(dir_link)

        restore_previous_session()

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "VDNA lab"))
        self.label_main_title.setText(_translate("MainWindow", "Virtual DNA lab"))
        self.label_write_browse_source_info.setText(_translate("MainWindow",
                                                               "Acceptable file formats ( .bngl ), Source BNGL file must include all parameters"))
        self.pushButton_write_browse_source.setText(_translate("MainWindow", "Browse"))
        self.label_write_browse_source.setText(_translate("MainWindow", "Source file"))
        self.label_write_list_created.setText(_translate("MainWindow", "List of ssDNA to be created"))
        self.pushButton_write_edit.setText(_translate("MainWindow", "Edit"))
        self.pushButton_write_delete.setText(_translate("MainWindow", "Delete"))
        self.pushButton_write_reset_all.setText(_translate("MainWindow", "Reset all"))
        self.checkBox_custom_file_name.setText(_translate("MainWindow", "Custom file name"))
        self.pushButton_write_import_sequence.setText(_translate("MainWindow", "Import from Species file"))
        self.pushButton_write_submit_without.setText(
            _translate("MainWindow", "Submit without existing ssDNAs sequence(s)"))
        self.label_write_allowed_char.setText(_translate("MainWindow", "Allowed characters (A, T, C, G)"))
        self.pushButton_write_submit_with.setText(
            _translate("MainWindow", "    Submit with existing ssDNAs sequence(s)   "))
        self.label_write_amount.setText(_translate("MainWindow", "Amount"))
        self.label_write_allowed_num.setText(_translate("MainWindow", "Numericals only"))
        self.pushButton_write_add_sequence.setText(_translate("MainWindow", "Add"))
        self.label_write_sequence.setText(_translate("MainWindow", "Desired ssDNA sequence"))
        self.pushButton_write_advanced.setText(_translate("MainWindow", "Advanced"))
        self.main_tab.setTabText(self.main_tab.indexOf(self.tab_write_bngl), _translate("MainWindow", "Write BNGL"))
        self.label_run_browse_source_info.setText(_translate("MainWindow",
                                                             "Acceptable file formats ( .bngl ), Source BNGL file must include all parameters"))
        self.pushButton_run_browse_source.setText(_translate("MainWindow", "Browse"))
        self.label_run_browse_source.setText(_translate("MainWindow", "Source file"))
        self.pushButton_run_browse_des.setText(_translate("MainWindow", "Browse"))
        self.label_run_browse_des_info.setText(
            _translate("MainWindow", "Dump file(s) willl be saved to this directory"))
        self.label_run_browse_des.setText(_translate("MainWindow", "Save simulation outputs directory"))
        self.pushButton_run_advanced.setText(_translate("MainWindow", "Advanced"))
        self.label_run_start_time.setText(_translate("MainWindow", "Start time"))
        self.label_run_n_dumps.setText(_translate("MainWindow", "# of dumps"))
        self.label_run_sim_end.setText(_translate("MainWindow", "Sim. end"))
        self.pushButton_run_run.setText(_translate("MainWindow", "Run"))
        self.main_tab.setTabText(self.main_tab.indexOf(self.tab_run_bngl), _translate("MainWindow", "Run BNGL"))
        self.label_visual_browse_source_info.setText(
            _translate("MainWindow", "Acceptable file formats ( .species, .0 )"))
        self.pushButton_visual_browse_source.setText(_translate("MainWindow", "Browse"))
        self.label_visual_browse_source.setText(_translate("MainWindow", "Source file"))
        self.pushButton_visual_view_source.setText(_translate("MainWindow", "View source species"))
        self.pushButton_visual_delete.setText(_translate("MainWindow", "Delete"))
        self.pushButton_visual_view_bngl.setText(_translate("MainWindow", "View result species"))
        self.pushButton_visual_view_comp.setText(_translate("MainWindow", "View complexes"))
        self.label_visual_history.setText(_translate("MainWindow", "History"))
        self.pushButton_visual_browse_des.setText(_translate("MainWindow", "Browse"))
        self.label_visual_browse_des_info.setText(
            _translate("MainWindow", "Analysis output files willl be saved to this directory"))
        self.label_visual_browse_des.setText(_translate("MainWindow", "Save analysis outputs directory"))
        self.pushButton_visual_advanced.setText(_translate("MainWindow", "Advanced"))
        self.pushButton_visual_visualize.setText(_translate("MainWindow", "Run visualization"))
        self.main_tab.setTabText(self.main_tab.indexOf(self.tab_visualize), _translate("MainWindow", "Visualize"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
