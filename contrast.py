from python_tools.correlator import Correlator
import click

@click.command()
@click.option('--handle', type=str, required=True)
@click.option('--data_filename', type=str, required=True)
@click.option('--box_size', type=float, required=True)
@click.option('--corr_type', type=str, required=True)
@click.option('--dim1_min', type=float, required=True)
@click.option('--dim1_max', type=float, required=True)
@click.option('--dim1_nbin', type=int, required=True)
@click.option('--dim2_min', type=float)
@click.option('--dim2_max', type=float)
@click.option('--dim2_nbin', type=int)
@click.option('--data_filename_2', type=str)

def run_contrast(
                 handle,
                 data_filename,
                 data_filename_2,
                 box_size,
                 corr_type,
                 dim1_min,
                 dim1_max,
                 dim1_nbin,
                 dim2_min,
                 dim2_max,
                 dim2_nbin
                 ):

    Correlator(handle=handle,
               data_filename=data_filename,
               data_filename_2=data_filename_2,
               box_size=box_size,
               corr_type=corr_type,
               dim1_min=dim1_min,
               dim1_max=dim1_max,
               dim1_nbin=dim1_nbin,
               dim2_min=dim2_min,
               dim2_max=dim2_max,
               dim2_nbin=dim2_nbin
               )


if __name__ == '__main__':
    run_contrast()