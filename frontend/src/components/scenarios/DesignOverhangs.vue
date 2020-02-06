<template lang="pug">
.page
  h1  {{ infos.title }}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Easy automated overhangs design for Golden Gate assembly",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')

  p.scenario-description.
    Find sets of compatible overhangs for Golden Gate assembly, to extend an
    assembly standard or decompose a sequence into compatible fragments

  //- learnmore Bla bla bla

  .form

    h4.formlabel What do you wish to do ?

    el-select(v-model='form.goal', placeholder='Select')
      el-option(label='Find a collection of compatible overhangs', value='overhangs_set')
      el-option(label='Decompose a sequence into compatible-ends fragments', value='sequence_decomposition')
      

    div(v-show='form.goal')

      div(v-show="form.goal === 'sequence_decomposition'")



        h4.formlabel Sequence

        collapsible(title='Examples')
          file-example(filename='phage_sequence.txt',
                       @input="function (e) { \
                         form.sequence = e; \
                         form.extremities = false; \
                         form.cutting_mode = 'equal'; \
                         form.n_fragments = 50; \
                        }",
                       fileHref='/static/file_examples/design_overhangs/phage_sequence.txt',
                       imgSrc='/static/file_examples/generic_logos/sequence_file.svg')
            p.
              A 50kb sequence of phage lambda.
              #[a(href="https://www.nature.com/articles/srep10655") Tsuge et al. ]
              managed to assemble this sequence in a one-step assembly of fifty 1kb
              fragments, digested by type-IIs enzymes. This required a very careful
              selection of enzyme overhangs.
          file-example(filename='sequence_with_specs.gb',
                         @input="function (e) { \
                           form.sequence = e; \
                           form.extremities = true; \
                           form.cutting_mode = 'features'; \
                          }",
                         fileHref='/static/file_examples/design_overhangs/sequence_with_specs.gb',
                         imgSrc='/static/file_examples/design_overhangs/sequence_with_specs.png')
              p.
                A 50kb sequence of phage lambda.
                #[a(href="https://www.nature.com/articles/srep10655") Tsuge et al. ]
                managed to assemble this sequence in a one-step assembly of fifty 1kb
                fragments, digested by type-IIs enzymes. This required a very careful
                selection of enzyme overhangs.

        files-uploader(v-model='form.sequence',
                       tip='Fasta or genbank. No file too large please :)',
                       :multiple='false')

        h4.formlabel
          span Sequence decomposition
          helper.
            If the provided sequence has features with a label "!cut",
            these features can be used as cutting zones. Otherwise the
            sequence will be cut into similar-length fragments')

        p
          el-radio(class='radio'
                   v-model='form.cutting_mode'
                   label='equal') Cut similar-length fragments

          el-radio(class='radio'
                   v-model='form.cutting_mode'
                   label='features') Cut in featured zones

        p
          el-checkbox(v-model='form.extremities') Consider overhangs on extremities
          helper.
            If this is checked, the overhangs at the suggested cut
            positions will also be compatible with the 4 basepairs
            at each extremity of the sequence

        p
          el-checkbox(v-model='form.allow_edits') Allow sequence edits
          helper.
            If this is checked, the sequence can be edited to make the
            decomposition possible. You can protect part of the sequence
            with features labeled "@DoNotModify" or "@EnforceTranslation".

        p(v-if="form.cutting_mode === 'equal'")
          span Number of fragments:
          el-input-number.inline(v-model="form.n_fragments",
                                 size="small",
                                 :min='2', :max='60')
        p Add sequences (e.g. enzyme sites) at the end of each fragment ? <br>
        p
          el-form(label-width="100px")
            el-form-item(label='Left flank:')
              el-input(v-model='left_flank_sequence', placeholder="e.g. CGTCTCA")
            el-form-item(label='Right flank:')
              el-input(v-model='right_flank_sequence', placeholder="e.g. TGAGACG")


      h4.formlabel Overhangs parameters

      p.inline(v-if="form.goal === 'overhangs_set'")
        span How many overhangs ?
        el-input-number.inline(v-show='!form.auto_overhangs'
                               v-model="form.n_overhangs", size="small",
                               :min='1', :max='50')
        el-checkbox.inline(v-model='form.auto_overhangs') as many as possible

      p.inline Differences between overhangs:
        el-input-number.inline(v-model="form.overhangs_differences", size="small",
                               :min='1', :max='4')
        //- helper(help='Number of basepairs by which every overhang pair should\
        //-              differ (this includes reverse complements).')
        helper.
          Number of basepairs by which every overhang pair should differ
          (this includes reverse complements).
      p.inline
        span Overhangs GC content
        el-slider.inline(range show-stops v-model='form.gc_content',
                         :max=100, :min=0, :step=25, :style="{width: '100px'}")
        span {{ form.gc_content[0] }}% to {{ form.gc_content[1] }}%:

      h4.formlabel {{form.specify_possible_overhangs ? 'Possible' : 'Forbidden'}} overhangs
      el-checkbox(v-model='form.specify_possible_overhangs') Specify possible overhangs rather than forbidden

      el-input(v-if='form.specify_possible_overhangs' type='textarea', :rows='4'
               v-model='form.possible_overhangs',
               placeholder='Enter possible overhangs, e.g. "ATTG, TTAC, ..."')
      el-input(v-else type='textarea', :rows='4'
               v-model='form.forbidden_overhangs',
               placeholder='Enter forbidden overhangs, e.g. "ATTG, TTAC, ..."')

      div
        h4.formlabel
          span External overhangs
          helper.
            The generated overhangs will be compatible with all the external overhangs
        el-input(type='textarea',
                 :rows='4',
                 v-model='form.mandatory_overhangs',
                 placeholder='Enter comma-separated overhangs, e.g. ATTG, TTAC, ...')

    backend-querier(:form='form',
                    :backendUrl='infos.backendUrl',
                    :validateForm='validateForm',
                    submitButtonText='Design',
                    v-show='form.goal',
                    v-model='queryStatus')

    p(v-if='queryStatus.polling.inProgress && queryStatus.polling.data.n_overhangs').center.
      Attempting to find {{queryStatus.polling.data.n_overhangs}} overhangs

    progress-bars(:bars='queryStatus.polling.data.bars', :order="['radius', 'interval']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')

    el-alert(v-show='queryStatus.requestError  && !queryStatus.polling.inProgress',
             :title="queryStatus.requestError",
             type="error",
             :closable="false")

  .results(v-if='!queryStatus.polling.inProgress && queryStatus.polling.data')

    p.center(v-if='!queryStatus.result.success').
      No solution found 😢. Maybe your parameters are too restrictive ?

    div(v-if='selected_overhangs')
      p We found a collection of {{ selected_overhangs.length }} overhangs:
      p.selected-overhangs
        span(v-for='overhang in selected_overhangs', :key='overhang') {{overhang + ', '}}

    download-button(v-if='queryStatus.result.zip_file',
                    :filedata='queryStatus.result.zip_file')

  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
// import filesuploader from '../../components/widgets/FilesUploader'
import tools from '../../tools'

var infos = {
  title: 'Design Golden Gate Overhangs',
  navbarTitle: 'Design Golden Gate Overhangs',
  path: 'design-overhangs',
  description: '',
  backendUrl: 'start/design_overhangs',
  icon: require('../../assets/images/overhangs.svg'),
  poweredby: ['goldenhinges', 'dnachisel']
}

export default {
  data () {
    return {
      overhangs: tools.cartesian(
        ['A', 'T', 'G', 'C'],
        ['A', 'T', 'G', 'C'],
        ['A', 'T', 'G', 'C'],
        ['A', 'T', 'G', 'C']
      ).map((l) => l.join('')),
      infos: infos,
      form: {
        goal: null,
        overhangs_differences: 1,
        gc_content: [25, 75],
        sequence: null,
        left_flank_sequence: '',
        right_flank_sequence: '',
        mandatory_overhangs: '',
        forbidden_overhangs: '',
        forbidden_overhangs_is_text: false,
        forbidden_overhangs_text: '',
        specify_possible_overhangs: false,
        possible_overhangs: '',
        cutting_mode: 'equal',
        n_overhangs: 10,
        n_fragments: 1,
        auto_overhangs: true,
        allow_edits: false,
        extremities: true
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  components: {
    learnmore
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if (!this.form.goal) {
        errors.push('Provide a goal !')
      } else if ((this.form.goal === 'sequence_decomposition') && (!this.form.sequence)) {
        errors.push('Provide a sequence !')
      }
      return errors
    }
  },
  computed: {
    selected_overhangs () {
      if (!this.queryStatus.result) {
        return null
      } else {
        return this.queryStatus.result.overhangs
      }
    },
    bars () {
      var data = this.queryStatus.polling.data
      if (!data) { return [] }
      return [
        {
          text: 'Cutting Interval',
          index: data.interval_ind,
          total: data.n_intervals
        }
      ]
    }
  }
}
</script>

<style lang='scss' scoped>

.form {
  margin: 50px auto;
  max-width: 500px;
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;
}

.el-checkbox {
  font-weight: normal;
  font-size: 16px !important;
}

.el-select {
  width: 100%
}

p.selected-overhangs {
  max-width: 90%;
  margin-left: 5%;

  span {
    font-family: 'Inconsolata', Courier;
    display: inline-block;
    margin: 0.2em 1em;
    padding: 3px;
    background-color: #eef3ff;
  }
}
</style>
