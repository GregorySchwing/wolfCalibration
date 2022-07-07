{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH -q gpu
#SBATCH --gres gpu:{{ gpus }}
{%- else %}
#SBATCH -q primary
#SBATCH --constraint=intel
{%- endif %}

#SBATCH -N 1
#SBATCH --mail-type=ALL

echo  "Running on host" hostname
echo  "Time is" date

conda activate nobel_gas

module load python/3.8
module swap gnu7 intel/2019

{% if gpus %}
module load cuda/11.0
#SBATCH --gres gpu:{{ gpus }}
{%- endif %}

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
