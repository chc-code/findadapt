import pickle
pw_script = os.path.dirname(os.path.realpath(__file__))
def get_all_miRNA_seq():
    fn = f'{pw_script}/miRNA.all.pre_miRNA_to_mature_miRNA.txt.gz'
    fn_pkl = f'{pw_script}/mirbase.all.pkl'
    # Accession	ID	Status	Sequence	Mature1_Acc	Mature1_ID	Mature1_Seq
    # MI0025295	bdi-MIR159b	UNCHANGED	UGGUUUUGAGGUGGAGCUCCUAUCAUUCCAAUGAAAGGUCUGUCAGAAGGGUGAUACAGCUGCUUGUUCAUGGUUCCCACUUAUCCUAUCUCCAUUAGAAAGGAGGAGGUAGGCUUGUGGUUUGCAUGACCGAGGAGCCGCAUCCACCCCUUGCCGACCGCUGUUUGGAUUGAAGGGAGCUCUGCAUCUUGAUCUUUC	MIMAT0030081	bdi-miR159b-5p.1	GAGCUCCUAUCAUUCCAAUGA

    # from column5, every 3 colmn is a group,  mature_acc, mature_id, mature_seq
    # each row may have different column number

    if os.path.exists(fn_pkl):
        with open(fn_pkl, 'rb') as f:
            return pickle.load(f)

    # usermap is for the user input converted to mature miRNA
    # idmap is the convert from miRBase ACC (MIxxx/MIMATxxx)  to the miRNA name/ID
    res = {'primary': {}, 'mature': {}, 'usermap': {}, 'idmap_pri': {}, 'idmap_mat': {}}

    # primary, key = primary name,  v = (pri_acc, pre_seq, [list of mature miRNA name])
    # mature, key = mature name, v = (mat_acc, mat_seq, primary_name, start_pos_at_pri_seq, end_pos_at_pri_seq)
    # idmap_pri, key = primary ACC, v = primary name
    # idmap_mat, key = mature_acc, v  mature name
    # usermap, convert/normalize the user input to get the most likely mature name

    def convert_seq(s):
        return s.upper().replace('U', 'T')

    with gzip.open(fn, 'rt') as f:
        f.readline()
        for i in f:
            i = i.strip()
            a = i.split('\t')
            pri_acc, pri_name, _, pri_seq = a[:4]
            pri_seq = convert_seq(pri_seq)
            pri_name = pri_name.lower()
            ires_pri = (pri_acc, pri_seq, [])
            res['primary'][pri_name] = ires_pri
            res['idmap_pri'][pri_acc] = pri_name

            mature_data = a[4:]
            for idx in range(int(len(mature_data)/3)):
                mat_acc, mat_name, mat_seq = mature_data[idx*3: idx*3 + 3]
                mat_name = mat_name.lower()
                mat_seq = convert_seq(mat_seq)
                start_pos_at_pri_seq = pri_seq.find(mat_seq)
                end_pos_at_pri_seq = start_pos_at_pri_seq + len(mat_seq)

                ires_mature = (mat_acc, mat_seq, pri_name, start_pos_at_pri_seq, end_pos_at_pri_seq)

                ires_pri[-1].append(mat_name)

                res['mature'][mat_name] = ires_mature
                res['idmap_mat'][mat_acc] = mat_name

                name_lite = re.sub('\-[35]p(\.\d+)?$', '', mat_name)
                name_lite2 = re.sub(r'[a-zA-Z]+$', '', name_lite)

                if name_lite != mat_name:
                    res['usermap'].setdefault(name_lite, []).append(mat_name)
                if name_lite2 != name_lite:
                    res['usermap'].setdefault(name_lite2, []).append(mat_name)

    with open(fn_pkl, 'wb') as o:
        pickle.dump(res, o)
    return res


def get_padding_seq(data, mat_data):
    res = {}
    pri_data_all = data['primary']
    for i in mat_data:
        mat_seq, pri_name, s, e = i[1:]
        pri_seq = pri_data_all[pri_name][1]
        seq_front = pri_seq[:s]
        seq_after = pri_seq[e:]
        res[mat_seq] = [seq_front, seq_after]
    return res


def prepare_refseq():
    # valid_organism = {'pma', 'pbv', 'mdm', 'zma', 'dvi', 'tgu', 'stu', 'oha', 'cgr', 'lja', 'ggo', 'ami', 'dre', 'aly', 'prd', 'ptc', 'cpi', 'cli', 'efu', 'ath', 'chi', 'cel', 'aca', 'cfa', 'ssc', 'dme', 'cja', 'hpo', 'ssa', 'gmo', 'bdi', 'cin', 'ppc', 'bmo', 'ocu', 'ptr', 'tca', 'pab', 'dno', 'oan', 'ppy', 'pal', 'cpo', 'eca', 'oni', 'osa', 'mtr', 'gma', 'rno', 'mml', 'bta', 'mdo', 'gga', 'mmu', 'hsa'}
    
    mature_data_all = get_all_miRNA_seq()
    valid_organism = {_.split('-', 1)[0] for _ in mir['mature'] if not _.startswith('mimat')} # 271 species

    for prefix in valid_organism:
        mat_data = [v for k, v in mature_data_all['mature'].items() if k[:3] == prefix]
        ref_data = get_padding_seq(mature_data_all, mat_data)
        with open(f'{pw_script}/{prefix}.refseq.pkl', 'wb') as o:
            pickle.dump(ref_data, o)
