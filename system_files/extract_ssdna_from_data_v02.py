

def extract_ssdna(source_path, c_from, dump_data):
    # reading raw data and return list with split data as needed
    # e.g N(b~A,5!3,3!4,W).N(b~A,5!25,3!26,W) 1
    # convert to...
    # ['N(b~A,5!3,3!4,W)', 'N(b~A,5!25,3!26,W)']
    # ['1']
    def data():
        def sep_data(f):
            li, ctl, syn = [], [], []
            for line in f:
                n_line = line.rstrip('\n').rstrip()
                if line.startswith('N'):
                    d = str(n_line).rsplit('  ', 1)[0]
                    ds = str(n_line).rsplit('  ', 1)[0].split('.')
                    ct = str(n_line).rsplit('  ', 1)[1]
                    li.append(ds)
                    ctl.append(int(ct))
                    syn.append(d)

            return [li, ctl, syn]

        if c_from == 'run_vis':
            with open(source_path, 'r') as f:
                file = f.readlines()
                f.close()

                get_return_lists = sep_data(file)
                return get_return_lists

        elif c_from == 'read_species':
            get_return_lists = sep_data(dump_data)
            return get_return_lists

    complexes_and_n = list(data())

    # converting raw data to ssDNA format
    # ['A', 'T'], if no compliments
    # ['A1', 'T1'], with compliments, letter and the compliment ID together
    def extract(d):
        b = ',5,'
        bm, lc = [], []
        for a in d:
            m, c = [], []
            bm.append(m)
            lc.append(c)
            for i in a:
                if b in i:
                    m.append(i)
                    c.append(d.index(a))

        w = []
        for en, an, in zip(bm, lc):
            f = []
            w.append(f)
            for e, a in zip(en, an):
                end = ',5!,'
                l, q = [], []
                f.append(l)
                q.append(e)
                while q[-1] != end:
                    for i in d[a]:
                        if q[-1] in i:
                            l_var = i.split(',')[2][2:]
                            letter = i.split(',')[0][4]
                            comp = i.split(',')[3][2:-1]
                            n = ',' + '5!' + l_var + ','
                            l.append(letter + comp)
                            q.append(n)
        return w

    data_set = extract(complexes_and_n[0])
    n_complexes = list(complexes_and_n[1])
    syns = list(complexes_and_n[2])

    if c_from == 'run_vis':
        return [data_set, n_complexes]

    elif c_from == 'read_species':
        ssdna_comp, ssdna_comp_n, ssdna_syn = [], [], []
        for comps, comp_n, syn in zip(data_set, n_complexes, syns):
            if comps not in ssdna_comp:
                ssdna_comp.append(comps)
                ssdna_comp_n.append(comp_n)
                ssdna_syn.append(syn)

            elif comps in ssdna_comp:
                comps_index = ssdna_comp.index(comps)
                ssdna_comp_n[comps_index] += comp_n

        ssdna_comp_ready = []
        for cc, cn, sy in zip(ssdna_comp, ssdna_comp_n, ssdna_syn):
            s_comp = []
            s_or_c = 2 if len(cc) > 1 else 2 \
                       if len(''.join([''.join([l for l in c if len(l) > 1]) for c in [e for e in cc]])) > 0 else 1
            l_only = [''.join([l[0] for l in c]) for c in [e for e in cc]]
            ss_sep = ' | '.join(l_only)
            s_comp.append(s_or_c)
            s_comp.append(cn)
            s_comp.append(ss_sep)
            s_comp.append(sy)
            ssdna_comp_ready.append(s_comp)

        return ssdna_comp_ready
