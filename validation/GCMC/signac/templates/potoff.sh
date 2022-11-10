{% extends "slurm.sh" %}

{% block header %}
{{- super () -}}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
{% set cpus = operations|map(attribute='directives.np')|sum %}

{% if gpus %}
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --exclude=ressrv7ai8111,ressrv8ai8111,ressrv13ai8111,ressrv14ai8111,ressrv15ai8111,res-lab34-ai8111,res-lab35-ai8111,res-lab41-ai8111,res-lab43-ai8111,reslab44ai8111
{%- else %}
#SBATCH -N 1
#SBATCH --exclude=reslab32ai8111,res-lab33-ai8111,res-lab34-ai8111
#SBATCH --mail-type=ALL
{%- endif %}


echo  "Running on host" $HOSTNAME
echo  "Time is" date

conda activate mosdef-study38

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
