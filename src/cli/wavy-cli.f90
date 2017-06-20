program wavycli

use mod_grid
use mod_domain
use mod_spectrum
use json_module
use flap,only:command_line_interface

implicit none

integer :: error
type(command_line_interface) :: cli

logical :: switch_value_domain
logical :: switch_value_grid
logical :: switch_value_spectrum

! initialize the CLI
call cli % init(progname='wavy-cli',version='0.1.0',authors='Wavebit Scientific LLC',&
                license='BSD-3-clause',&
                description='Command line interface to libwavy',&
                examples=["wavy-cli --help"])

! define sub-commands
call cli % add_group(group='new',description='create a new instance')
call cli % add_group(group='load',description='loads an existing instance from json')

! define command line arguments
call cli % add(group='new',switch='--spectrum',switch_ab='-s',&
               help='Create new spectrum',required=.false.,def='.false.',&
               act='store_true',error=error)

call cli % add(group='new',switch='--domain',switch_ab='-d',&
               help='Create new domain',required=.false.,def='.false.',&
               act='store_true',error=error)

call cli % add(group='new',switch='--grid',switch_ab='-g',&
               help='Create new grid',required=.false.,def='.false.',&
               act='store_true',error=error)

call cli % get(group='new',switch='--spectrum',val=switch_value_spectrum)
call cli % get(group='new',switch='--domain',val=switch_value_domain)
call cli % get(group='new',switch='--grid',val=switch_value_grid)

endprogram wavycli
