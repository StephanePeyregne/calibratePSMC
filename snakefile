import yaml

# Load configuration
with open("config.yaml") as f:
    config = yaml.safe_load(f)


individual = config["individual"]
data = config["data"]
simulation_code = config["simulation_code"].format(individual=individual)
genome_length = config["genome_length"].format(individual=individual)
filters = config["filters"].split()
scrm2psmc_script = config["scripts"]["scrm2psmc"]
psmc_binary = config["scripts"]["psmc"]



THETA = [x / 100.0 for x in range(100, 201, 5)]
RHO = [x / 100.0 for x in range(100, 201, 5)]


rule simulations:
    input:
        expand("simulations/{individual}_{theta}_{rho}.likelihood_round20", individual=individual, theta=THETA, rho=RHO)


rule generate_simulation:
    input:
        simulation_code
    params:
        u = "{theta}",
        r = "{rho}"
    output:
        "simulations/{individual}_{theta}_{rho}.sh"
    shell:
        """awk -v rho={params.r} -v theta={params.u} '{{printf "scrm 2 1 -SC abs -p 10 -t "theta*$5" -r "rho*$7; for (i=8; i<=NF; i++) printf " "$i}}' {input} > {output}"""


rule run_simulation:
    input:
        script = "simulations/{individual}_{theta}_{rho}.sh"
    output:
        "simulations/{individual}_{theta}_{rho}.scrm"
    shell:
        """bash {input.script} > {output}"""


rule scrm2psmc:
    input:
        scrm = "simulations/{individual}_{theta}_{rho}.scrm",
        length = genome_length,
        filters = filters,
        script = scrm2psmc_script
    output:
        "simulations/{individual}_{theta}_{rho}.psmcfa"
    shell:
        """python3 {input.script} -g {input.length} -m <(zcat {input.filters}) {input.scrm} > {output}"""


rule psmc_simulation:
    input:
        simulated_data = "simulations/{individual}_{theta}_{rho}.psmcfa",
        psmc = psmc_binary
    output:
        "simulations/{individual}_{theta}_{rho}.psmc"
    shell:
        """{input.psmc} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output} {input.simulated_data}"""


rule psmc2demography:
    input:
        "simulations/{individual}_{theta}_{rho}.psmc"
    output:
        "simulations/{individual}_{theta}_{rho}_round20.demography"
    shell:
        """grep "^PA" {input} | awk '{{if (NR==21) {{printf $2; for (i=3; i<=NF; i++) printf " "$i}}}}' > {output}"""


rule likelihood:
    input:
        data = data,
        demography = "simulations/{individual}_{theta}_{rho}_round20.demography",
        psmc = psmc_binary
    output:
        "simulations/{individual}_{theta}_{rho}.likelihood"
    shell:
        """{input.psmc} -N1 -i {input.demography} -o {output} {input.data}"""


rule likelihood_table:
    input:
        "simulations/{individual}_{theta}_{rho}.likelihood"
    params:
        individual = individual
    output:
        "simulations/{individual}_{theta}_{rho}.likelihood_round20"
    shell:
        """ cat {input} | awk -v u={wildcards.theta} -v r={wildcards.rho} -v ind={params.individual} 'BEGIN{{to_print=0}}{{if ($1=="RD" && $2==1) {{to_print=1}}; if ($1=="LK" && to_print==1) {{print ind","u","r","$2; to_print=0}} }}' > {output} """


def summary(wildcards):
    return expand("simulations/{individual}_{theta}_{rho}.likelihood_round20", individual=individual, theta=THETA, rho=RHO)

rule summary_table:
    input:
        summary
    output:
        "simulations/likelihood_table.csv"
    shell:
        """ cat {input} >> {output} """


