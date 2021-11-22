import subprocess

def vcf_to_list(vcf_f, increment = 1, start = 0):
    from vcf import Reader
    vcf_reader = Reader(open(vcf_f, 'r'))
    vcf_data = []
    count = 0
    next_count = start
    for record in vcf_reader:
        if count == next_count:
            vcf_data.append(record)
            next_count += increment
        count += 1
    vcf_data = sorted(sorted(vcf_data, key = lambda x: x.POS), key = lambda x: x.CHROM)
    samples = [x.sample for x in vcf_data[0].samples]
    return (vcf_data, samples)

def extract_var_from_bed(bed, bcf, out_dir):
    print("Extracting variants from " + bcf + " based on positions in " + bed)
    tmp_vcf = out_dir + "/variants"
    subprocess.call(("vcftools", "--bcf", bcf, "--recode", "--bed", bed, "--out", tmp_vcf))
    print("Successfuly extracted variants into " + tmp_vcf + ".recode.vcf")
    return tmp_vcf + ".recode.vcf"
