<template lang="pug">
.page
  h1  {{ infos.title }}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Easy automated overhangs design for Golden Gate assembly",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')

  p.scenario-description Find sets of compatible overhangs for your assembly problem.

  //- learnmore Bla bla bla

  .form
    h4.formlabel Enzyme Sites to choose from:
    .example(@click="form.sites_to_include = exampleEnzymes.join(', ')"
             style='cursor: pointer; color: grey; margin: 0.5em;') Example (click me)
    el-input(type='textarea', :rows='7'
              v-model='form.sites_to_include',
              placeholder='Enter enzyme names, e.g. "EcoRI, PvuI, ..."')
    p
      el-checkbox(v-model='form.enforce_unique_sites') A site should appear at most once
    h4.formlabel Sites you absolutely want in the sequence: 
    el-input(type='textarea', :rows='5'
              v-model='form.mandatory_sites',
              placeholder='Enter enzyme names, e.g. "EcoRI, PvuI, ..."')
    h4.formlabel Forbidden enzyme sites: 
    el-input(type='textarea', :rows='3'
              v-model='form.forbidden_sites',
              placeholder='Enter enzyme names, e.g. "EcoRI, PvuI, ..."')

    backend-querier(:form='form',
                    :backendUrl='infos.backendUrl',
                    :validateForm='validateForm',
                    submitButtonText='Design',
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

    p {{queryStatus.result.message}}
    download-button(:filedata='queryStatus.result.genbank_file',
                    text='Download Genbank')
    center
      iframe(:src="queryStatus.result.figure_data"
              name="stacked_arrays.pdf"
              style="width: 90%; max-width: 1000px; height:300px;"
              frameborder="0")

  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Design Sites Arrays',
  navbarTitle: 'Design Sites Arrays',
  path: 'design-sites-arrays',
  description: '',
  backendUrl: 'start/design_sites_arrays',
  icon: require('../../assets/images/design_sites_arrays.svg'),
  poweredby: ['zymp']
}

export default {
  data () {
    return {
      infos: infos,
      form: {
        sites_to_include: '',
        mandatory_sites: '',
        forbidden_sites: '',
        enforce_unique_sites: false
      },
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      exampleEnzymes: [
        'AccI', 'AclI', 'AflII', 'AflIII', 'AgeI', 'ApaLI', 'AseI',
        'AvaI', 'BamHI', 'BanII', 'BlnI', 'BmtI', 'BsmI', 'BssHII',
        'DdeI', 'DraI', 'Eco47III', 'EcoRI', 'EcoRV', 'HindII',
        'HindIII', 'HinfI', 'HpaI', 'KpnI', 'MfeI', 'MluI',
        'MspA1I', 'MunI', 'NaeI', 'NcoI', 'NdeI', 'NheI', 'NotI',
        'NsiI', 'NspI', 'PstI', 'PvuI', 'PvuII', 'SacI', 'SacII',
        'SalI', 'ScaI', 'SfaNI', 'SnaBI', 'SpeI', 'SphI', 'SspI',
        'StyI', 'VspI', 'XhoI', 'XmaI', 'ZraI'
      ]
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
      if (!this.form.sites_to_include.length) {
        errors.push('Provide at least some enzymes to choose from!')
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
