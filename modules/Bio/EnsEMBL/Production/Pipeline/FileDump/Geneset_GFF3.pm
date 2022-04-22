=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Production::Pipeline::FileDump::Geneset_GFF3;

use strict;
use warnings;
no  warnings 'redefine';
use base qw(Bio::EnsEMBL::Production::Pipeline::FileDump::Base_Filetype);

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::CDS;
use Bio::EnsEMBL::Utils::IO::GFFSerializer;
use Path::Tiny;
use Bio::EnsEMBL::MiscFeature;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    data_type            => 'genes',
    file_type            => 'gff3',
    per_chromosome       => 0,
    feature_types        => ['Gene', 'Transcript'],
    include_utr          => 1,
    header               => 1,
    gt_gff3_exe          => 'gt gff3',
    gt_gff3validator_exe => 'gt gff3validator',
  };
}

sub run {
  my ($self) = @_;

  my $data_type      = $self->param_required('data_type');
  my $per_chromosome = $self->param_required('per_chromosome');
  my $feature_types  = $self->param_required('feature_types');
  my $header         = $self->param_required('header');
  my $filenames      = $self->param_required('filenames');

  my $dba = $self->dba;

  my %adaptors;
  foreach my $feature_type (@$feature_types) {
    $adaptors{$feature_type} = $dba->get_adaptor($feature_type);
  }

  my $mca = $dba->get_adaptor('MetaContainer');
  my $providers = $mca->list_value_by_key('assembly_provider.name') || ();
  my $provider = join("; ", @$providers);

  my ($chr, $non_chr, $non_ref) = $self->get_slices($dba);

  my $filename = $$filenames{$data_type};

  if ($per_chromosome && scalar(@$chr)) {
    $self->print_to_file($chr, 'chr', $filename, '>', $header, $dba, \%adaptors, $provider);
    if (scalar(@$non_chr)) {
      $self->print_to_file($non_chr, 'non_chr', $filename, '>>', $header, $dba, \%adaptors, $provider);
    }
  } else {
    $self->print_to_file([@$chr, @$non_chr], undef, $filename, '>', $header, $dba, \%adaptors, $provider);
  }

  if (scalar(@$non_ref)) {
    my $non_ref_filename = $self->generate_non_ref_filename($filename);
    path($filename)->copy($non_ref_filename);
    $self->print_to_file($non_ref, undef, $non_ref_filename, '>>', $header, $dba, \%adaptors, $provider);
  }

  my $output_filenames = $self->param('output_filenames');
  foreach (@$output_filenames) {
    $self->tidy($_);
    $self->validate($_);
  }
}

sub print_to_file {
  my ($self, $slices, $region, $filename, $mode, $header, $dba, $adaptors, $provider) = @_;

  $header = $mode eq '>' ? 1 : 0 unless defined $header;
  my $serializer = $self->gff3_serializer($filename, $mode, $header, $slices, $dba);

  my $non_chr_serializer;
  if ($region && $region eq 'non_chr') {
    my $non_chr_filename = $self->generate_non_chr_filename($filename);
    $non_chr_serializer = $self->gff3_serializer($non_chr_filename, $mode, 1, $slices, $dba);
  }

  while (my $slice = shift @{$slices}) {
    my $chr_serializer;
    if ($region && $region eq 'chr') {
      my $chr_filename = $self->generate_chr_filename($filename, $slice);
      $chr_serializer = $self->gff3_serializer($chr_filename, $mode, 1, $slices, $dba);
    }

    $slice->source($provider) if defined $provider;
    $serializer->print_feature($slice);
    if ($region && $region eq 'chr') {
      $chr_serializer->print_feature($slice);
    } elsif ($region && $region eq 'non_chr') {
      $non_chr_serializer->print_feature($slice);
    }

    foreach my $feature_type (keys %$adaptors) {
      my $features = $self->fetch_features($feature_type, $$adaptors{$feature_type}, $slice);
      $serializer->print_feature_list($features);
      if ($region && $region eq 'chr') {
        $chr_serializer->print_feature_list($features);
      } elsif ($region && $region eq 'non_chr') {
        $non_chr_serializer->print_feature_list($features);
      }
    }
  }
}

sub gff3_serializer {
  my ($self, $filename, $mode, $header, $slices, $dba) = @_;

  my $fh = path($filename)->filehandle($mode);
  my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($fh);
  if ($header) {
    $serializer->print_main_header($slices, $dba);
  }

  return $serializer;
}

sub fetch_features {
  my ($self, $feature_type, $adaptor, $slice) = @_;

  my @features = @{$adaptor->fetch_all_by_Slice($slice)};
  $adaptor->clear_cache();
  if ($feature_type eq 'Transcript') {
    my $exon_features = $self->exon_features(\@features);
    push @features, @$exon_features;
  }

  return \@features;
}

sub exon_features {
  my ($self, $transcripts) = @_;
  my $include_utr = $self->param_required('include_utr');

  my @cds_features;
  my @exon_features;
  my @utr_features;
  my @codon_features;

  foreach my $transcript (@$transcripts) {
    push @cds_features, @{ $transcript->get_all_CDS(); };
    push @exon_features, @{ $transcript->get_all_ExonTranscripts() };
    if ($include_utr) {
      push @utr_features, @{ $transcript->get_all_five_prime_UTRs()};
      push @utr_features, @{ $transcript->get_all_three_prime_UTRs()};
    }

    push @codon_features, get_start_stop_codon_features($transcript);
    push @codon_features, make_seleno_features($transcript);
  }

  return [@exon_features, @cds_features, @utr_features, @codon_features];
}

sub tidy {
  my ($self, $filename) = @_;
  my $gt_gff3_exe = $self->param_required('gt_gff3_exe');

  $self->assert_executable($gt_gff3_exe);

  my $temp_filename = "$filename.sorted";

  my $cmd = "$gt_gff3_exe -tidy -sort -retainids -fixregionboundaries -force -o $temp_filename $filename";
  my ($rc, $output) = $self->run_cmd($cmd);

  if ($rc) {
    my $msg =
      "GFF3 tidying failed for '$filename'\n".
      "Output: $output";
    unlink $temp_filename;
    $self->throw($msg);
  } else {
    path($temp_filename)->move($filename);
  }
}

sub validate {
  my ($self, $filename) = @_;
  my $gt_gff3validator_exe = $self->param_required('gt_gff3validator_exe');

  $self->assert_executable($gt_gff3validator_exe);

  my $cmd = "$gt_gff3validator_exe $filename";
  my ($rc, $output) = $self->run_cmd($cmd);

  if ($rc) {
    my $msg =
      "GFF3 file '$filename' failed validation with $gt_gff3validator_exe\n".
      "Output: $output";
    $self->throw($msg);
  }
}

sub get_start_stop_codon_features {
  my ($transcript) = @_;

  my @codon_features;

  my @start_codons = make_start_codon_features($transcript);
  my @stop_codons = make_stop_codon_features($transcript);

  my $translation = $transcript->translation();
  return (()) unless (defined($translation));

  my ($has_start, $has_end) = check_start_and_stop($transcript);

  foreach my $exon (@{$transcript->get_all_Exons}) {
    if ($translation && $exon == $translation->start_Exon && $has_start) {
      push(@codon_features, @start_codons);
    }
    if ($translation && $exon == $translation->end_Exon && $has_end) {
      push(@codon_features, @stop_codons);
    }
  }

  return @codon_features;
}

sub check_start_and_stop {
  my ($transcript) = @_;

  return (0, 0) unless (defined($transcript->translation));
  my ($has_start, $has_end);

  my @attrib = @{$transcript->get_all_Attributes('cds_start_NF')};
  $has_start = (scalar(@attrib) == 1 && $attrib[0]->value() == 1) ? 0 : 1;
  @attrib = @{$transcript->get_all_Attributes('cds_end_NF')};
  $has_end = (scalar(@attrib) == 1 && $attrib[0]->value() == 1) ? 0 : 1;

  return (0, 0) unless ($has_start || $has_end);

  my $cds_seq = uc($transcript->translateable_seq);
  my $start_seq = substr($cds_seq, 0, 3);

  my @exons = @{$transcript->get_all_translateable_Exons};
  my $last_exon = $exons[$#exons];
  my $phase = $last_exon->end_phase;
  $phase = 0 if ($phase == -1);
  my $end_seq = substr($cds_seq, -(3-$phase));

  my ($attrib) = @{$transcript->slice()->get_all_Attributes('codon_table')};
  my $codon_table_id = $attrib->value() if (defined($attrib));
  $codon_table_id ||= 1;
  my $codon_table = Bio::Tools::CodonTable->new(-id => $codon_table_id);

  $has_start = 0 unless ($codon_table->is_start_codon($start_seq));
  $has_end = 0 unless ($codon_table->is_ter_codon($end_seq));

  return ($has_start, $has_end);
}

sub make_start_codon_features {
  my ($transcript) = @_;

  return (()) if (!$transcript->translation);

  my @pepgencoords = $transcript->pep2genomic(1,1);

  return (()) if (scalar(@pepgencoords) > 2);
  return (()) unless ($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
  return (()) unless ($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate'));

  my @translateable_exons = @{$transcript->get_all_translateable_Exons};
  my @start_codon_features;

  foreach my $pepgencoord (@pepgencoords) {
    my $misc_feature = Bio::EnsEMBL::MiscFeature->new(
      -seqname => $transcript->seq_region_name,
      -start => $pepgencoord->start,
      -end   => $pepgencoord->end,
      -strand => $translateable_exons[0]->strand,
    );
    my $attributes = {
      'so_term' => 'start_codon',
      'parent_so_term' => 'transcript',
      'version' => $transcript->version,
      'parent' => $transcript->stable_id,
      'source' => $transcript->source
    };
    while (my ($key, $value) = each %{$attributes}) {
      my $attr = Bio::EnsEMBL::Attribute->new(
        -CODE => $key,
        -NAME => $key,
        -VALUE => $value
      );
      $misc_feature->add_Attribute($attr);
    }
    push @start_codon_features, $misc_feature;
  }

  if ($translateable_exons[0]->strand == 1) {
    @start_codon_features = sort {$a->start <=> $b->start } @start_codon_features;
  } else {
    @start_codon_features = sort {$b->start <=> $a->start } @start_codon_features;
  }

  return @start_codon_features;
}

sub make_stop_codon_features {
  my ($transcript) = @_;

  return (()) if (!$transcript->translation);

  my $cdna_endpos = $transcript->cdna_coding_end;
  my @pepgencoords = $transcript->cdna2genomic($cdna_endpos-2,$cdna_endpos);

  return (()) if(scalar(@pepgencoords) > 2);
  return (()) unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate'));
  return (()) unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate'));

  my @translateable_exons = @{$transcript->get_all_translateable_Exons};
  my @stop_codon_features;

  foreach my $pepgencoord (@pepgencoords) {
    my $misc_feature = Bio::EnsEMBL::MiscFeature->new(
      -seqname => $transcript->seq_region_name,
      -start => $pepgencoord->start,
      -end   => $pepgencoord->end,
      -strand => $translateable_exons[0]->strand,
    );
    my $attributes = {
      'so_term' => 'stop_codon',
      'parent_so_term' => 'transcript',
      'version' => $transcript->version,
      'parent' => $transcript->stable_id,
      'source' => $transcript->source
    };
    while (my ($key, $value) = each %{$attributes}) {
      my $attr = Bio::EnsEMBL::Attribute->new(
        -CODE => $key,
        -NAME => $key,
        -VALUE => $value
      );
      $misc_feature->add_Attribute($attr);
    }
    push @stop_codon_features, $misc_feature;
  }

  if ($translateable_exons[0]->strand == 1) {
    @stop_codon_features = sort {$a->start <=> $b->start } @stop_codon_features;
  } else {
    @stop_codon_features = sort {$b->start <=> $a->start } @stop_codon_features;
  }

  return @stop_codon_features;
}

sub make_seleno_features {
  my ($transcript) = @_;

  my @seleno_codon_features;

  my $selenos = check_selenos($transcript);
  if (exists($selenos->{$transcript->stable_id})) {
    my ($intrans, $instop) = (0, 0);
    my $slice_offset = $transcript->slice->start-1;
    my @translateable_exons = @{$transcript->get_all_translateable_Exons} if $transcript->translation;

    foreach my $exon (@{$transcript->get_all_Exons}) {
      if ($transcript->translation && $exon == $transcript->translation->start_Exon){
        $intrans = 1;
      }

      if ($intrans) {
        my $cds_exon = shift @translateable_exons;
        return(()) if (!$cds_exon);

        my $exon_start = $cds_exon->start;
        my $exon_end = $cds_exon->end;

        if (!$instop){
          my $cds_start = ($exon_start + $slice_offset);
          my $cds_end = ($exon_end + $slice_offset);

          foreach my $seleno_mod (split(",", $selenos->{$transcript->stable_id})) {
            push(@seleno_codon_features, create_seleno_attributes($transcript, $seleno_mod, $cds_start, $cds_end));
          }
        }
      }
    }
  }

  return @seleno_codon_features;
}

sub check_selenos {
  my ($transcript) = @_;

  my %selenos;
  my @selenos_attributes = ();

  my $translation = $transcript->translation();
  if ($translation) {
    my $attributes = $translation->get_all_Attributes('_selenocysteine');

    foreach my $attr (@{$attributes}) {
      $attr->value =~ /^(\d+).+/;
      if ($1) {
        push(@selenos_attributes, $1);
      }
    }
  }

  if (scalar(@selenos_attributes)) {
    $selenos{$transcript->stable_id} = join(",", sort {$a<=>$b} @selenos_attributes);
  }

  return \%selenos;
}

sub create_seleno_attributes {
  my ($transcript, $aa_start_pos, $cds_start, $cds_end) = @_;

  my $transcript_mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);
  my @coords = $transcript_mapper->pep2genomic($aa_start_pos, $aa_start_pos);
  my @seleno_codon_features;

  foreach my $coord (@coords) {
    if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      if ($coord->start >= $cds_start && $coord->end <= $cds_end) {
        my $misc_feature = Bio::EnsEMBL::MiscFeature->new(
          -seqname => $transcript->seq_region_name,
          -start => $coord->start,
          -end   => $coord->end,
          -strand => $transcript->strand,
        );
        my $attributes = {
          'so_term' => 'stop_codon_redefined_as_selenocysteine',
          'parent_so_term' => 'transcript',
          'version' => $transcript->version,
          'parent' => $transcript->stable_id,
          'source' => $transcript->source
        };
        while (my ($key, $value) = each %{$attributes}) {
          my $attr = Bio::EnsEMBL::Attribute->new(
            -CODE => $key,
            -NAME => $key,
            -VALUE => $value
          );
          $misc_feature->add_Attribute($attr);
        }
        push @seleno_codon_features, $misc_feature;
      }
    }
  }

  return @seleno_codon_features;
}

# The serializer will put everything from 'summary_as_hash' into
# the GFF3 attributes. But we don't want everything, so redefine
# here, to explicitly specify our attributes.

sub Bio::EnsEMBL::Gene::summary_as_hash {
  my $self = shift;
  my %summary;

  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $self->source;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->display_id;
  $summary{'gene_id'}         = $self->display_id;
  $summary{'version'}         = $self->version;
  $summary{'Name'}            = $self->external_name;
  $summary{'biotype'}         = $self->biotype;
  $summary{'description'}     = $self->description;
  
  return \%summary;
}

sub Bio::EnsEMBL::Transcript::summary_as_hash {
  my $self = shift;
  my %summary;

  my $parent_gene = $self->get_Gene();

  $summary{'seq_region_name'}          = $self->seq_region_name;
  $summary{'source'}                   = $self->source;
  $summary{'start'}                    = $self->seq_region_start;
  $summary{'end'}                      = $self->seq_region_end;
  $summary{'strand'}                   = $self->strand;
  $summary{'id'}                       = $self->display_id;
  $summary{'Parent'}                   = $parent_gene->stable_id;
  $summary{'transcript_id'}            = $self->display_id;
  $summary{'version'}                  = $self->version;
  $summary{'Name'}                     = $self->external_name;
  $summary{'biotype'}                  = $self->biotype;
  my $havana_transcript = $self->havana_transcript();
  $summary{'havana_transcript'}        = $havana_transcript->display_id() if $havana_transcript;
  $summary{'havana_version'}           = $havana_transcript->version() if $havana_transcript;
  $summary{'ccdsid'}                   = $self->ccds->display_id if $self->ccds;
  $summary{'transcript_support_level'} = $self->tsl if $self->tsl;
  $summary{'tag'}                      = 'basic' if $self->gencode_basic();

  # Check for seq-edits
  my $seq_edits = $self->get_all_SeqEdits();
  my @seq_edits;
  foreach my $seq_edit (@$seq_edits) {
    my ($start, $end, $alt_seq) =
      ($seq_edit->start, $seq_edit->end, $seq_edit->alt_seq);

    my $note = "This transcript's sequence has been ";
    if ($alt_seq eq '') {
      $note .= "annotated with a deletion between positions $start and $end";
    } elsif ($end < $start) {
      $note .= "annotated with an insertion before position $start: $alt_seq";
    } else {
      $note .= "replaced between positions $start and $end: $alt_seq";
    }

    push @seq_edits, $note;
  }
  $summary{'Note'} = \@seq_edits if scalar(@seq_edits);

  return \%summary;
}

sub Bio::EnsEMBL::Exon::summary_as_hash {
  my $self = shift;
  my %summary;

  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'exon_id'}         = $self->display_id;
  $summary{'version'}         = $self->version;
  $summary{'constitutive'}    = $self->is_constitutive;

  return \%summary;
}

sub Bio::EnsEMBL::CDS::summary_as_hash {
  my $self = shift;
  my %summary;

  $summary{'seq_region_name'} = $self->seq_region_name;
  $summary{'source'}          = $self->transcript->source;
  $summary{'start'}           = $self->seq_region_start;
  $summary{'end'}             = $self->seq_region_end;
  $summary{'strand'}          = $self->strand;
  $summary{'phase'}           = $self->phase;
  $summary{'id'}              = $self->translation_id;
  $summary{'Parent'}          = $self->transcript->display_id;
  $summary{'protein_id'}      = $self->translation_id;
  $summary{'version'}         = $self->transcript->translation->version;

  return \%summary;
}

sub Bio::EnsEMBL::MiscFeature::summary_as_hash {
  my $self = shift;
  my %summary;

  $summary{'seq_region_name'} = $self->seqname;
  $summary{'source'}          = $self->get_scalar_attribute('source');
  $summary{'start'}           = $self->start;
  $summary{'end'}             = $self->end;
  $summary{'strand'}          = $self->strand;
  $summary{'id'}              = $self->get_scalar_attribute('parent');
  $summary{'Parent'}          = $self->get_scalar_attribute('parent');
  $summary{'version'}         = $self->get_scalar_attribute('version');
  $summary{'so_term'}         = $self->get_scalar_attribute('so_term');
  $summary{'parent_so_term'}  = $self->get_scalar_attribute('parent_so_term');

  return \%summary;
}

1;
