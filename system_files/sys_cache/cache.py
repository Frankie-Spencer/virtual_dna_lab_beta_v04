import os

write_b_source = 'wbs'
run_b_source, run_b_des = 'rbs', 'rbd'
visual_b_source, visual_b_des = 'vbs', 'vbd'
import_seq = 'isf'

what_loc = {write_b_source: 'system_files/sys_cache/_write_source_loc.cache',
            run_b_source: 'system_files/sys_cache/_run_source_loc.cache',
            run_b_des: 'system_files/sys_cache/_run_des_loc.cache',
            import_seq: 'system_files/sys_cache/_write_import_loc.cache',
            visual_b_source: 'system_files/sys_cache/_visual_source_loc.cache',
            visual_b_des: 'system_files/sys_cache/_visual_des_loc.cache'}


def read_s_loc(what):
    with open(what_loc[what], 'r') as cf:
        li = cf.readline()
        cf.close()
        return li


def write_s_loc(what, l):
    with open(what_loc[what], 'w') as cf:
        cf.writelines(l)
        cf.close()


def read_d_loc(what):
    with open(what_loc[what], 'r') as cf:
        li = cf.readline()
        cf.close()
        return li


def write_d_loc(what, l):
    with open(what_loc[what], 'w') as cf:
        cf.writelines(l)
        cf.close()


def write_browser_loc(l):
    with open('system_files/sys_cache/_browser_loc.cache', 'w') as cf:
        cf.writelines(l)
        cf.close()


def read_browser_loc():
    with open('system_files/sys_cache/_browser_loc.cache', 'r') as cf:
        li = cf.readline()
        cf.close()
        return li
