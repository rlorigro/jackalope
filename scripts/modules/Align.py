from multiprocessing import Pool
import subprocess
import sys
import os


def run_minigraph(output_directory, gfa_path, fasta_path, n_threads, args_override=None):
    output_path = os.path.join(output_directory, "alignments.gaf")

    # minigraph \
    # -cx lr \
    # -o reads_vs_graph.gaf \
    # graph.gfa \
    # reads.fasta \
    # args = ["minigraph", "-c", "-x", "lr", "-o", output_path, gfa_path, fasta_path]

    if args_override is None:
        args = [
            "minigraph",
            "-c",
            "-g", str(10000),
            "-k", str(14),
            "-f", "0.25",
            "-r", "1000,20000",
            "-n", "3,3",
            "-p", str(0.5),
            # "-j", str(0.85),  # <-- this alone causes horrific slowdown in some regions, no idea why
            "-x", "lr",
            "-t", str(n_threads),
            "-o", output_path,
            gfa_path,
            fasta_path]
    else:
        args = \
            ["minigraph"] + \
            args_override + \
            ["-o", output_path,
             gfa_path,
             fasta_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(args)+'\n')

        try:
            p1 = subprocess.run(args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    sys.stderr.write("done\n")

    return output_path


def run_panaligner(output_directory, gfa_path, fasta_path, n_threads, args_override=None):
    output_path = os.path.join(output_directory, "reads_vs_graph.gaf")

    # minigraph \
    # -cx lr \
    # -o reads_vs_graph.gaf \
    # graph.gfa \
    # reads.fasta \
    # args = ["minigraph", "-c", "-x", "lr", "-o", output_path, gfa_path, fasta_path]

    if args_override is None:
        args = [
            "PanAligner",
            "-c",
            "-g", str(10000),
            "-k", str(14),
            "-f", "0.25",
            "-r", "1000,20000",
            "-n", "3,3",
            "-p", str(0.5),
            # "-j", str(0.85),  # <-- this alone causes horrific slowdown in some regions, no idea why
            "-x", "lr",
            "-t", str(n_threads),
            "-o", output_path,
            gfa_path,
            fasta_path]
    else:
        args = \
            ["PanAligner"] + \
            args_override + \
            ["-o", output_path,
             gfa_path,
             fasta_path]

    with open(output_path, 'a') as file:
        sys.stderr.write(" ".join(args)+'\n')

        try:
            p1 = subprocess.run(args, stdout=file, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            sys.stderr.write("Status : FAIL " + '\n' + (e.stderr.decode("utf8") if e.stderr is not None else "") + '\n')
            sys.stderr.flush()
            return False
        except Exception as e:
            sys.stderr.write(str(e))
            return False

    sys.stderr.write("done\n")

    return output_path


