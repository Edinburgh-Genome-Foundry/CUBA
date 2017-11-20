<template lang="pug">
.page
  h1  {{ infos.title }}

  img.icon.center-block(slot='title-img', :src='infos.icon')

  p.center Find sets of compatible overhangs for your assembly problem.

  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Design Golden Gate overhangs online:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")

  learnmore Bla bla bla

  .form

    h4.formlabel What do you wish to design ?

    el-select(v-model='form.goal', placeholder='Select')
      el-option(v-for='item in goal_options', :label='item.label',
                :value='item.value', :key='item.value')

    div(v-show='form.goal')
      div(v-show="form.goal === 'sequence_decomposition'")


        h4.formlabel Sequence

        examples-dialog

        filesuploader(v-model='form.sequence',
                      text="Drop a single file (or click to select)",
                      help='Fasta or genbank. No file too large please :)',
                      :multiple='false')

        h4.formlabel
          span Sequence decomposition
          helper(help='If the provided sequence has features with a label "!cut",\
                       these features can be used as cutting zones. Otherwise the\
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
          helper(help="If this is checked, the overhangs at the suggested cut\
                       positions will also be compatible with the 4 basepairs\
                       at each extremity of the sequence")

        p
          el-checkbox(v-model='form.allow_edits') Allow sequence edits
          helper(help='If this is checked, the sequence can be edited to make the\
                       decomposition possible. You can protect part of the sequence\
                       with features labeled "@DoNotModify" or "@EnforceTranslation".')

        p(v-if="form.cutting_mode === 'equal'")
          span Number of fragments:
          el-input-number.inline(v-model="form.n_fragments",
                                 size="small",
                                 :min=2, :max=60)

      h4.formlabel Overhangs parameters

      p.inline(v-if="form.goal === 'overhangs_set'")
        span How many overhangs ?
        el-input-number.inline(v-show='!form.auto_overhangs'
                               v-model="form.n_overhangs", size="small",
                               :min=1, :max=50)
        el-checkbox.inline(v-model='form.auto_overhangs') as many as possible

      p.inline Differences between overhangs:
        el-input-number.inline(v-model="form.overhangs_differences", size="small",
                               :min=1, :max=4)
        helper(help='Number of basepairs by which every overhang pair should\
                     differ (this includes reverse complements).')

      p.inline
        span Overhangs GC content
        el-slider.inline(range show-stops v-model='form.gc_content',
                         :max=100, :min=0, :step=25, :style="{width: '100px'}")
        span {{ form.gc_content[0] }}% to {{ form.gc_content[1] }}%:

      h4 Forbidden overhangs

      el-input(type='textarea',
               :rows='4'
               v-model='form.forbidden_overhangs',
               placeholder='Enter comma-separated overhangs, e.g. ATTG, TTAC, ...')

      div(v-show="form.goal === 'overhangs_set'")
        h4 External overhangs
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

    progress-bars(:bars='bars')

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
import filesuploader from '../../components/widgets/FilesUploader'
import tools from '../../tools'

var infos = {
  title: 'Design Golden Gate Overhangs',
  navbarTitle: 'Design Golden Gate Overhangs',
  path: 'design-overhangs',
  description: '',
  backendUrl: 'start/design_overhangs',
  icon: require('assets/images/overhangs.svg'),
  poweredby: ['goldenhinges', 'dnachisel']
}

export default {
  data () {
    return {
      overhangs: tools.cartesian(['A', 'T', 'G', 'C'],
                                 ['A', 'T', 'G', 'C'],
                                 ['A', 'T', 'G', 'C'],
                                 ['A', 'T', 'G', 'C']).map((l) => l.join('')),
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
      },
      goal_options: [
        {
          label: 'A collection of compatible overhangs',
          value: 'overhangs_set'
        },
        {
          label: 'A sequence decomposition, with compatible overhangs',
          value: 'sequence_decomposition'
        }
      ]
    }
  },
  components: {
    filesuploader,
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

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
}

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

.el-checkbox.inline {
  margin-left: 15px;
}


.el-slider.inline {
  display: inline-block;
  margin-left: 20px;
  margin-right: 20px;
  margin-bottom: -12px;
}

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
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
