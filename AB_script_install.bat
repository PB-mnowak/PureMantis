start https://www.microsoft.com/en-us/p/python-39/9p7qfqmjrfp7
pause

pip install -r requirements.txt

mkdir "%HOMEDRIVE%%HOMEPATH%%\Desktop\Sequencing"
copy "P:\_research group folders\AB Antibodies and Phage Display\_CDR analysis script\cdr_analysis_script_pball.bat" "%HOMEDRIVE%%HOMEPATH%%\Desktop\Sequencing\"
mkdir "%HOMEDRIVE%%HOMEPATH%%\Desktop\Sequencing\scf_files"
pause