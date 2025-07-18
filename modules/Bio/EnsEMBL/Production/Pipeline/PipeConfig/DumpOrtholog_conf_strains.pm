=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::Production::Pipeline::PipeConfig::DumpOrtholog_conf_strains;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::PipeConfig::DumpOrtholog_conf');

use Bio::EnsEMBL::Hive::Version;

sub default_options {
	my ($self) = @_;

	return {
		%{ $self->SUPER::default_options() },

		species_config => {
      '1' => {
              compara     => 'multi',
              source      => 'mus_musculus',
              antispecies => ['mus_musculus','mus_spretus','mus_pahari','mus_caroli'],
              taxons      => ['Mus'],
              division    => 'EnsemblVertebrates',
              homology_types =>
                ['ortholog_one2one', 'apparent_ortholog_one2one'],
      },
	  }
	};
}

1;
